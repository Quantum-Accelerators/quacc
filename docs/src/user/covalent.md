# Using Covalent

Here, we will show how to use [Covalent](https://github.com/AgnostiqHQ/covalent) to construct, dispatch, and monitor workflows in Quacc. If you prefer to use Jobflow and/or FireWorks, skip to the ["Using Quacc with Jobflow and FireWorks"](jobflow.md) section of the documentation.

## Pre-Requisites

Make sure you completed the "Covalent Setup" section of the documentation. If you haven't done so already, run `covalent start` prior to starting these examples to start the Covalent server and UI.

Additionally, you should read the Covalent [First Experiment](https://covalent.readthedocs.io/en/latest/getting_started/first_experiment/index.html) guide to get a sense of how Covalent works. Namely, you should understand the [Covalent Basics](https://covalent.readthedocs.io/en/latest/concepts/basics.html) of the `Electron` and `Lattice` objects, which describe individual compute tasks and workflows, respectively.

In Covalent, the `@ct.lattice` decorator indicates that the function is a workflow, and the `@ct.electron` decorator indicates that the function is a job (i.e. an individual compute task). If you plan to use a job scheduling system like Slurm, you can think of each `Electron` as an individual Slurm job. If `Electron` objects are wrapped by a `Lattice`, they will only be executed once the job is dispatched. Conversely, if you did not include the `@ct.lattice` decorator, all the `Electron` objects would behave as normal Python functions.

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
    result1 = relax_job(atoms)

    # Define Job 2
    result2 = static_job(result1["atoms"])

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

Here, you can see that it is quite trivial to set up a workflow using the recipes within Quacc. We define the full workflow as a `Lattice` object and stitch together the individual workflow steps. By default, all Quacc recipe functions are `Electron` objects, so we didn't need to use the `@ct.electron` decorator around the functions here.

Covalent will also automatically construct a directed acyclic graph of the inputs and outputs for each calculation to determine which jobs are dependent on one another and the order the jobs should be run. In this example, Covalent will know not to run `job2` until `job1` has completed. While that may seem obvious here, this can be quite powerful when you start to construct more complex workflows and start to dispatch individual `Electron` calculations on remote compute resources.

### Running a Simple Parallel Workflow

Now let's consider a similar but nonetheless distinct example. Here, we will define a workflow where we will carry out two EMT structure relaxations, but the two jobs are not dependent on one another. In this example, Covalent will know that it can run the two jobs separately, and even if Job 1 were to fail, Job 2 would still progress.

```python
import covalent as ct
from ase.build import bulk, molecule
from quacc.recipes.emt.core import relax_job

# Define workflow
@ct.lattice
def workflow(atoms1, atoms2):

    # Define two separate relaxation jobs
    result1 = relax_job(atoms1)
    result2 = relax_job(atoms2)

    return [result1, result2]

# Define two Atoms objects
atoms1 = bulk("Cu")
atoms2 = molecule("N2")

# Dispatch the workflow to the Covalent server
dispatch_id = ct.dispatch(workflow)(atoms1, atoms1)

# Fetch the results from the server
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

### Running Workflows with Complex Connectivity

Now that you understand how to write and run simple workflows consisting of Quacc recipes, let's consider a more complex example.

For this example, let's consider a toy scenario where we wish to relax a bulk Cu structure, carve all possible slabs, and then run a new relaxation calculation on each slab.

In Quacc, there are two types of recipes: 1) individual compute tasks that are functions; 2) defined workflows that are classes. Here, we are interested in importing a workflow, so it will be instantiated slightly differently from the prior examples. See the example below:

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import BulkToSlabsFlow

@ct.lattice
def workflow(atoms):
    relaxed_bulk = relax_job(atoms)
    relaxed_slabs = BulkToSlabsFlow(static_electron=None).run(relaxed_bulk["atoms"])
    return relaxed_slabs

atoms = bulk("Cu")
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

Here, we have imported to the `BulkToSlabsFlow` class, which is instantiated with optional parameters and is applied to an `Atoms` object. Here, for demonstration purposes, we specify the `static_electron=None` option to only do a relaxation after the slabs are made. All we have to do to use the workflow is wrap it inside a `@ct.lattice` decorator. We don't need to use a `@ct.electron` decorator here because the `BulkToSlabsFlow` class is already calling `Electron` objects internally. As a general rule, all classes in Quacc are workflows that can be transformed into a `Lattice`, and all functions are compute tasks that can be transformed into `Electron` objects.

If you want to understand what is going on underneath the hood, it is worth checking out the source code. Because the number of slabs is not pre-determined, this recipe is using a Covalent feature called a [Sublattice](https://covalent.readthedocs.io/en/latest/concepts/basics.html?highlight=sublattice#sublattice) that enables [dynamic workflows](https://covalent.readthedocs.io/en/latest/developer/patterns/dynamic_workflow.html?highlight=dynamic%20workflow). Of course, if you don't plan to develop new dynamic workflows, you don't actually need to know this. You can just import the recipe that makes use of this feature and away you go.

## Setting Executors

By defualt, Covalent will run all `Electron` tasks on your local machine using the [`DaskExecutor`](https://covalent.readthedocs.io/en/latest/api/executors/dask.html). This is a parameter that you can control. For instance, you may want to define the executor to be based on [Slurm](https://covalent.readthedocs.io/en/latest/api/executors/slurm.html) to submit a job to an HPC cluster. The example below highlights how one can change the executor for individual jobs.

If you want to use the same executor for all the `Electron` objects in a `Lattice`, you can pass the `executor` keyword argument to the `@ct.lattice` decorator, as shown below.

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

@ct.lattice(executor="local")
def workflow(atoms):

    result1 = relax_job(atoms)
    result2 = static_job(result1["atoms"])

    return result2

atoms = bulk("Cu")
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

If you construct workflow from scratch, you can also set the executor for individual `Electron` objects by passing the `executor` keyword argument to the `@ct.electron` decorator. In Quacc, the `Electron` executor options can be modified after they are imported as follows:

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

@ct.lattice
def workflow(atoms):
    job1 = relax_job
    job1.electron_object.metadata["executor"] = "dask"

    job2 = static_job
    job2.electron_object.metadata["executor"] = "local"

    output1 = job1(atoms)
    output2 = job2(output1["atoms"])
    return output2

atoms = bulk("Cu")
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

## Learn More

That ends the Covalent section of the documentation. If you want to learn more about Covalent, you can read the [Covalent Documentation](https://covalent.readthedocs.io/en/latest/index.html).