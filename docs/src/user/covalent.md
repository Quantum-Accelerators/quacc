# Using Quacc with Covalent

Here, we will show how to use [Covalent](https://github.com/AgnostiqHQ/covalent) to construct, dispatch, and monitor workflows in Quacc.

```{note}
If you prefer to use a workflow engine other than Covalent, then refer to the ["Alternate Workflow Engines"](advanced/alt_workflows/index.md) section of the documentation.
```

## Pre-Requisites

Make sure you completed the ["Covalent Setup"](../install/covalent.md) section of the documentation. Additionally, you should learn about the main [Covalent Concepts](https://docs.covalent.xyz/docs/user-documentation/concepts/concepts-index), namely the [`Electron`](https://docs.covalent.xyz/docs/user-documentation/concepts/covalent-basics#electron) and [`Lattice`](https://docs.covalent.xyz/docs/user-documentation/concepts/covalent-basics#lattice) objects, which describe individual compute tasks and workflows, respectively.

In Covalent, the `@ct.lattice` decorator indicates that the function is a workflow, and the `@ct.electron` decorator indicates that the function is a job (i.e. an individual compute task). If you plan to use a job scheduling system like Slurm, you can think of each `Electron` as an individual Slurm job. For some minimal working examples of how to write your own Covalent workflows and how they compare to other workflow tools, refer to the [Worfklow Engine Comparison Guide](advanced/alt_workflows/comparison.md).

All `Electron` and `Lattice` objects behave as normal Python functions when the necessary arguments are supplied. However, if the `ct.dispatch` command is used, the workflow will be dispatched to the Covalent server for execution and monitoring.

## Examples

```{hint}
If you haven't done so yet, make sure you started the Covalent server with `covalent start` in the command-line.
```

### Running a Simple Serial Workflow

We will now try running a simple workflow where we relax a bulk Cu structure using EMT and take the output of that calculation as the input to a follow-up static calculation with EMT.

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

# Define the workflow
@ct.lattice
def workflow(atoms):

    # Define Job 1
    result1 = relax_job(atoms)

    # Define Job 2, which takes the output of Job 1 as input
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

You can see that it is quite trivial to set up a workflow using the recipes within Quacc. We define the full workflow as a `Lattice` object that stitches together the individual workflow steps. The {obj}`.emt.core.relax_job` and {obj}`.emt.core.static_job` were both already defined with a `@ct.electron` decorator, so they will be interpreted by Covalent as `Electron` objects.

Covalent will also automatically construct a directed acyclic graph of the inputs and outputs for each calculation to determine which jobs are dependent on one another and the order the jobs should be run. In this example, Covalent will know not to run `job2` until `job1` has completed successfully.

The job will be dispatched to the Covalent server with the [`ct.dispatch`](https://docs.covalent.xyz/docs/user-documentation/concepts/covalent-basics#dispatch) command, which takes in the workflow function and the input arguments to the workflow. The [`ct.get_result`](https://docs.covalent.xyz/docs/user-documentation/concepts/covalent-basics#result) command is used to fetch the results from the server.

```{note}
Because the workflow is only sent to the server with `ct.dispatch`, calling `workflow(atoms)` would run the workflow as if Covalent were not being used at all.
```

![Covalent UI](../_static/user/tutorial1.jpg)

### Running a Simple Parallel Workflow

Now let's consider a similar but nonetheless distinct example. Here, we will define a workflow where we will carry out two EMT structure relaxations, but the two jobs are not dependent on one another. In this example, Covalent will know that it can run the two jobs separately, and even if Job 1 were to fail, Job 2 would still progress.

```python
import covalent as ct
from ase.build import bulk, molecule
from quacc.recipes.emt.core import relax_job

# Define workflow
@ct.lattice
def workflow(atoms1, atoms2):

    # Define two independent relaxation jobs
    result1 = relax_job(atoms1)
    result2 = relax_job(atoms2)

    return {"result1": result1, "result2": result2}

# Define two Atoms objects
atoms1 = bulk("Cu")
atoms2 = molecule("N2")

# Dispatch the workflow to the Covalent server
dispatch_id = ct.dispatch(workflow)(atoms1, atoms2)

# Fetch the results from the server
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

![Covalent UI](../_static/user/tutorial2.jpg)

### Running Workflows with Complex Connectivity

For this example, let's consider a toy scenario where we wish to relax a bulk Cu structure, carve all possible slabs, and then run a new relaxation calculation on each slab (with no static calculation at the end). This is an example of a dynamic workflow.

In Quacc, there are two types of recipes: individual compute tasks with the suffix `_job` and pre-made multi-step workflows with the suffix `_flow`. Here, we are interested in importing a pre-made workflow. Refer to the example below:

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow

@ct.lattice
def workflow(atoms):
    relaxed_bulk = relax_job(atoms)
    relaxed_slabs = bulk_to_slabs_flow(relaxed_bulk["atoms"], slab_static_electron=None)

    return relaxed_slabs

atoms = bulk("Cu")
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

We have imported the {obj}`.emt.slabs.bulk_to_slabs_flow` function, which takes an `Atoms` object along with several optional parameters. For demonstration purposes, we specify the `slab_static_electron=None` option to do a relaxation but disable the static calculation on each slab. All we have to do to define the workflow is wrap it inside a `@ct.lattice` decorator.

Due to the dynamic nature of `bulk_to_slabs_flow`, the number of returned slabs will be dependent on the input `Atoms` object. The pattern for creating a dynamic workflow in Covalent is called a ["sublattice"](https://docs.covalent.xyz/docs/user-documentation/concepts/covalent-arch/covalent-sdk#sublattice). The sublattice, which is really just a fancy name for a sub-workflow within a larger workflow, and its individual compute tasks can also be viewed in the Covalent UI.

```{hint}
You don't need to set `wait=True` in practice. Once you call `ct.dispatch`, the workflow will begin running. The `ct.get_result` function is used to fetch the workflow status and results from the server.
```

![Covalent UI](../_static/user/tutorial3.gif)

## Setting Executors

By default, Covalent will run all `Electron` tasks on your local machine using the [`DaskExecutor`](https://docs.covalent.xyz/docs/user-documentation/api-reference/executors/dask). This is a parameter that you can control. For instance, you may want to define the executor to be based on [Slurm](https://docs.covalent.xyz/docs/user-documentation/api-reference/executors/slurm) to submit a job to an HPC cluster. The example below highlights how one can change the executor.

### Setting Executors via the Lattice Object

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

### Setting Executors via the Electron Objects

The individual `Electron` executor options can be modified after they are imported as follows:

```python
import covalent as ct
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

@ct.lattice
def workflow(atoms):
    job1 = relax_job
    job1.electron_object.executor = "dask"

    job2 = static_job
    job2.electron_object.executor = "local"

    output1 = job1(atoms)
    output2 = job2(output1["atoms"])
    return output2

atoms = bulk("Cu")
dispatch_id = ct.dispatch(workflow)(atoms)
result = ct.get_result(dispatch_id, wait=True)
print(result)
```

```{hint}
If you are defining your own workflow functions to use, you can also set the executor for individual `Electron` objects by passing the `executor` keyword argument to the `@ct.electron` decorator.
```

## Configuring Executors

Refer to the [executor documentation](https://docs.covalent.xyz/docs/features/executor-plugins/exe) for instructions on how to configure Covalent for your desired machines.

By default, the `workdir` for the `Dask` (default) and `local` executors is set to `~/.cache/covalent/workdir`. This is where any files generated at runtime will be stored. You can change both of these parameters to the directories of your choosing by editing the Covalent configuration file directly or via the `ct.set_config()` command.

For submitting jobs to [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/) from your local machine, an example `SlurmExecutor` configuration with support for an [`sshproxy`](https://docs.nersc.gov/connect/mfa/#sshproxy)-based multi-factor authentication certificate might look like the following:

```python
n_nodes = 1
n_cores_per_node = 48

executor = ct.executor.SlurmExecutor(
    username="YourUserName",
    address="perlmutter-p1.nersc.gov",
    ssh_key_file="~/.ssh/nersc",
    cert_file="~/.ssh/nersc-cert.pub",
    remote_workdir="$SCRATCH",
    conda_env="quacc",
    options={
        f"nodes": {n_nodes},
        "qos": "debug",
        "constraint": "cpu",
        "account": "YourAccountName",
        "job-name": "quacc",
        "time": "00:10:00",
    },
    prerun_commands=[
        "export COVALENT_CONFIG_DIR=$SCRATCH",
        f"export QUACC_VASP_PARALLEL_CMD='srun -N {n_nodes} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores'",
    ],
    use_srun=False,
)
```

```{important}
The `SlurmExecutor` must have `use_srun=False` in order for ASE-based calculators to be launched appropriately.
```

## Learn More

That ends the Covalent section of the documentation. If you want to learn more about Covalent, you can read the [Covalent Documentation](https://docs.covalent.xyz/docs/). Please refer to the Covalent [Discussion Board](https://github.com/AgnostiqHQ/covalent/discussions) for any Covalent-specific questions.
