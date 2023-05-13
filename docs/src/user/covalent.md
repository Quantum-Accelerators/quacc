# Using Covalent

Here, we will show how to use [Covalent](https://github.com/AgnostiqHQ/covalent) to construct, dispatch, and monitor workflows in QuAcc.

## Pre-Requisites

Make sure you completed the ["Covalent Setup"](covalent.md) section of the documentation. If you haven't done so already, run `covalent start` prior to starting these examples.

Additionally, you should read the Covalent [First Experiment](https://covalent.readthedocs.io/en/latest/getting_started/first_experiment/index.html) guide to get a sense of how Covalent works. Namely, you should understand the [Covalent Basics](https://covalent.readthedocs.io/en/latest/concepts/basics.html) of the `Electron` and `Lattice` objects, which describe individual compute tasks and workflows, respectively.

## Running a Simple Serial Workflow 

Here, we will try running a simple workflow where we relax a bulk Cu structure using EMT and take the output of that calculation as the input to a follow-up static calculation with EMT. This example is shown below.

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

# Define the workflow
@ct.lattice
def workflow(atoms):

    # Define Job 1
    job1 = ct.electron(static_job)

    # Define Job 2
    job2 = ct.electron(relax_job)

    # Run Job 1
    result1 = job1(atoms)

    # Run Job 2 on the output Atoms from Job 1
    result2 = job2(result1["atoms"])

    return result2

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Dispatch the workflow to the Covalent server
# with the bulk Cu Atoms object as the input
dispatch_id = ct.dispatch(workflow)(atoms)

# Fetch the result from the server
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

Here, you can see that it is quite trivial to set up a workflow using the recipes within QuAcc. In Covalent, the `@ct.lattice` decorator indicates that the function is a workflow, and the `@ct.electron` decorator indicates that the function is a job (i.e. an individual compute task). Here, we have simply used the shorthand of `ct.electron(<function>)`.

If `Electron` objects are wrapped by a `Lattice` function, they will only be executed once the job is dispatched. Conversely, if you did not include the `@ct.lattice` decorator, all the `Electron` objects would behave as normal Python functions.

Covalent will also automatically construct a directed acyclic graph of the inputs and outputs for each calculation to determine which jobs are dependent on one another and the order the jobs should be run. In this example, Covalent will know not to run `job2` until `job1` has completed. This can be quite powerful when you start to construct more complex workflows.

### Running a Simple Parallel Workflow

Now let's consider a similar but nonetheless distinct example. Here, we will define a workflow where we will carry out two EMT structure relaxations, but the two jobs are not dependent on one another. In this example, Covalent will know that it can run the two jobs separately of one another, and even if Job 1 were to fail, Job 2 would still progress.

```python
import covalent as ct
from ase.build import bulk, molecule
from quacc.recipes.emt.core import relax_job

# Define workflow
@ct.lattice
def workflow(atoms1, atoms2):

    # Define and Run two separate relaxation jobs
    result1 = ct.electron(relax_job)(atoms1)
    result2 = ct.electron(relax_job)(atoms2)

    return [result1, result2]

# Define two Atoms objects
bulk = bulk("Cu")
molecule = molecule("N2")

# Dispatch the workflow to the Covalent server
dispatch_id = ct.dispatch(workflow)(bulk,molecule)

# Fetch the results from the server
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

Note that in this example, we have used `ct.electron(<function>)(*args, **kwargs)` in a single line for brevity, but you can always split it into two like the prior example if you prefer.

### Running Workflows with Complex Connectivity

Now that you understand how to write and run simple workflows consisting of QuAcc recipes, let's consider a more complex example.

For this example, let's consider a toy scenario where we wish to relax a bulk Cu structure, carve all possible slabs, and then run a new relaxation calculation on each slab. Because the number of slabs is not pre-determined, we will need to use a Covalent feature called a [Sublattice](https://covalent.readthedocs.io/en/latest/concepts/basics.html?highlight=sublattice#sublattice), which is illustrated below.


Since the topic of dynamic workflows is a bit complex, QuAcc provides a pre-packaged workflow for this exact scenario. This is shown below.

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.vasp.slabs import BulkToSlabsFlow

@ct.lattice
def workflow(atoms)
    return BulkToSlabsFlow(static_electron=None).run(atoms)

atoms = bulk("Cu")
dispatch_id = ct.dispatch(workflow)(atoms)
results = ct.get_result(dispatch_id, wait=True)
print(results)
```

Now you know how to run dynamic workflows from QuAcc and a little bit about what is happening underneath the hood.

## Setting Executors

By defualt, Covalent will run all `Electron` tasks on your local machine using the [`DaskExecutor`](https://covalent.readthedocs.io/en/latest/api/executors/dask.html). This is a parameter that you can control. For instance, you will often want to define the executor to be based on [Slurm](https://covalent.readthedocs.io/en/latest/api/executors/slurm.html) to submit a job to an HPC cluster. The example below highlights how one can change the executor for individual jobs.

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import static_job

@ct.electron(executor="dask")
def relax_electron(atoms):
    return relax_job(atoms)

@ct.electron(executor="local")
def static_electron(atoms):
    return static_job(atoms)

@ct.lattice
def workflow(atoms):
    output1 = relax_electron(atoms)
    output2 = static_electron(output1["atoms"])
    return output2

atoms = bulk("Cu")
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

Here, we have chosen to use the `@ct.electron` decorator around a function instead of the `ct.electron(<function>)` shorthand so that we could specify the `executor` keyword argument.

Alternatively, this can be done more compactly by directly setting the `metadata` attribute of the `Electron` object, as shown below.

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

@ct.lattice
def workflow(atoms):
    job1 = ct.electron(relax_job)
    job1.electron_object.metadata["executor"] = "dask"

    job2 = ct.electron(static_job)
    job2.electron_object.metadata["executor"] = "local"

    output1 = job1(atoms)
    output2 = job2(output1["atoms"])
    return output2

atoms = bulk("Cu")
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

If you want to use the same executor for all the `Electron` objects in a `Lattice`, you can pass the `executor` keyword argument to the `@ct.lattice` decorator instead of the individual `Electron` objects.

## Learn More

That ends the Covalent section of the documentation. If you want to learn more about Covalent, you can read the [Covalent Documentation](https://covalent.readthedocs.io/en/latest/index.html).