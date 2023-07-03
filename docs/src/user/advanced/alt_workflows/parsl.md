# Using Quacc with Parsl

## Introduction

[Parsl](https://github.com/Parsl/parsl) is a Python program developed at Argonne National Laboratory, the University of Chicago, and the University of Illinois to easily write parallel workflows that can be dispatched on distributed compute resources. Like Jobflow+FireWorks, it can be used in place of Covalent, if preferred.

## Pre-Requisites

Make sure you completed the ["Parsl Setup"](../../../install/advanced/alt_workflows/parsl.md) section of the installation instructions. Additionally, you should read the Parsl documentation's ["Quick Start"](https://parsl.readthedocs.io/en/stable/quickstart.html) to get a sense of how Parsl works. Namely, you should understand the concept of a `@python_app` and `@join_app`, which describe individual compute tasks and dynamic job tasks, respectively.

```{seealso}
For a more detailed tutorial on how to use Parsl, refer to the ["Parsl Tutorial"](https://parsl.readthedocs.io/en/stable/1-parsl-introduction.html) and the even more detailed ["Parsl User Guide"](https://parsl.readthedocs.io/en/stable/userguide/index.html).
```

## Examples

```{hint}
If you haven't loaded your Parsl config, you must do that first so Parsl can construct the job dependency graph. For testing purposes, you simply can run `import parsl` followed by `parsl.load()` before starting the examples below, which will enable jobs to run on your local machine.
```

### Running a Simple Serial Workflow

We will first try running a simple workflow where we relax a bulk Cu structure using EMT and take the output of that calculation as the input to a follow-up static calculation with EMT. Note that all dependencies need to be defined within the `@python_app` definition.

```python
from parsl import python_app
from ase.build import bulk

# Define the Python apps
@python_app
def relax_app(atoms):
    from quacc.recipes.emt.core import relax_job

    return relax_job(atoms)

@python_app
def static_app(atoms):
    from quacc.recipes.emt.core import static_job

    return static_job(atoms)

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Call App 1
future1 = relax_app(atoms)

# Call App 2, which takes the output of App 1 as input
future2 = static_app(future1)

# Print result
print(future2.result())
```

You can see that it is quite trivial to set up a Parsl workflow using the recipes within Quacc. We define the full workflow as a function that stitches together the individual `@python_app` workflow steps.

```{note}
The use of `.result()` serves to block any further calculations from running until it is resolved. Calling `.result()` also returns the function output as opposed to the `AppFuture` object.
```

### Running a Simple Parallel Workflow

Now let's consider a similar but nonetheless distinct example. Here, we will define a workflow where we will carry out two EMT structure relaxations, but the two jobs are not dependent on one another. In this example, Parsl will know that it can run the two jobs in parallel, and even if Job 1 were to fail, Job 2 would still progress.

```python
from parsl import python_app
from ase.build import bulk, molecule

# Define the Python app
@python_app
def relax_app(atoms):
    from quacc.recipes.emt.core import relax_job

    return relax_job(atoms)

# Define two Atoms objects
atoms1 = bulk("Cu")
atoms2 = molecule("N2")

# Define two independent relaxation jobs
future1 = relax_app(atoms1)
future2 = relax_app(atoms2)

# Print the results
print(future1.result(), future2.result())
```

If you monitor the output, you'll notice that the two jobs are being run in parallel, whereas they would be run sequentially if you were not using Parsl.

### Running Workflows with Complex Connectivity

#### The Inefficient Way

For this example, let's consider a toy scenario where we wish to relax a bulk Cu structure, carve all possible slabs, and then run a new relaxation calculation on each slab (with no static calculation at the end). This is an example of a dynamic workflow.

In Quacc, there are two types of recipes: individual compute tasks with the suffix `_job` and pre-made multi-step workflows with the suffix `_flow`. Here, we are interested in importing a pre-made workflow. Refer to the example below:

```python
from parsl import python_app
from ase.build import bulk

@python_app
def relax_app(atoms):
    from quacc.recipes.emt.core import relax_job

    return relax_job(atoms)

@python_app
def bulk_to_slabs_app(atoms):
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    return bulk_to_slabs_flow(atoms, slab_static_electron=None)

# Define the Atoms object
atoms = bulk("Cu")

# Define the workflow
future1 = relax_app(atoms)
future2 = bulk_to_slabs_app(future1)

# Print the results
print(future2.result())
```

When running a Covalent-based workflow like {obj}`.emt.slabs.bulk_to_slabs_flow` above, the entire function will run as a single compute task even though it is composed of several individual sub-tasks. If these sub-tasks are compute-intensive, this might not be the most efficient use of resources.

#### The Efficient Way

Quacc fully supports Parsl-based workflows to resolve this limitation. For example, the workflow above can be equivalently run as follows using the Parsl-specific {obj}`.emt.parsl.slabs.bulk_to_slabs_app` workflow:

```python
from parsl import python_app
from ase.build import bulk
from quacc.recipes.emt.parsl.slabs import bulk_to_slabs_flow

# Define the Python App
@python_app
def relax_app(atoms):
    from quacc.recipes.emt.core import relax_job

    return relax_job(atoms)

# Define the Atoms object
atoms = bulk("Cu")

# Define the workflow
future1 = relax_app(atoms)
future2 = bulk_to_slabs_flow(future1, slab_static_app=None)

# Print the results
print(future2.result())
```

In this example, all the individual tasks and sub-tasks are run as separate jobs, which is more efficient. By comparing {obj}`.emt.parsl.slabs.bulk_to_slabs_app` with its Covalent counterpart {obj}`.emt.slabs.bulk_to_slabs_flow`, you can see that the two are extremely similar such that it is often straightforward to [interconvert](comparison.md) between the two.

```{note}
We didn't need to wrap `bulk_to_slabs_flow` with a decorator because it returns an `AppFuture`. This is also why we call `.result()` on it.
```

## Job Management

Out-of-the-box, Parsl will run on your local machine. However, in practice you will probably want to run your Parsl workflows on HPC machines.

```{note}
If you are just starting out, try running some test calculations locally first. Then come back and set up the relevant configuration files for your desired machines.
```

### Configuring Executors

To configure Parsl for the high-performance computing environment of your choice, refer to the executor [Configuration](https://parsl.readthedocs.io/en/stable/userguide/configuring.html) page in the Parsl documentation.

For [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/), example `HighThroughputExecutor` configurations can be found in the [NERSC Documentation](https://docs.nersc.gov/jobs/workflow/parsl/). A simple one is reproduced below that allows for job submission from the login node. This example will create a single Slurm job that will run one `PythonApp` at a time on a single node and is good for testing out some of the examples above.

```python
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SimpleLauncher
from parsl.providers import SlurmProvider

config = Config(
    max_idletime=120,
    executors=[
        HighThroughputExecutor(
            label="quacc_HTEX",
            max_workers=1,
            provider=SlurmProvider(
                account="MyAccountName",
                nodes_per_block=1,
                scheduler_options="#SBATCH -q debug -C cpu",
                worker_init="source activate quacc",
                walltime="00:10:00",
                cmd_timeout=120,
                launcher = SimpleLauncher(),
            ),
        )
    ],
)
```

The individual arguments are as follows:

- `max_idletime`: The maximum amount of time (in seconds) to allow the executor to be idle before the Slurm job is cancelled.
- `label`: A label for the executor instance, used during file I/O.
- `max_workers`: Maximum number of workers to allow on a node.
- `SlurmProvider()`: The provider to use for job submission. This can be changed to `LocalProvider()` if you wish to have the Parsl process run on a compute node rather than the login node.
- `account`: Your NERSC account name.
- `nodes_per_block`: The number of nodes to request per job. By default, all cores on the node will be requested (seetting `cores_per_node` will override this).
- `scheduler_options`: Any additional `#SBATCH` options can be included here.
- `worker_init`: Commands to run before the job starts, typically used for activating a given Python environment.
- `walltime`: The maximum amount of time to allow the job to run in `HH:MM:SS` format.
- `cmd_timeout`: The maximum time to wait (in seconds) for the job scheduler info to be retrieved/sent.
- `launcher`: The type of Launcher to use. Note that `SimpleLauncher()` must be used instead of the commonly used `SrunLauncher()` to allow Quacc subprocesses to launch their own `srun` commands.

Unlike some other workflow engines, Parsl (by default) is built for "jobpacking" where the allocated nodes continually pull in new workers (until the walltime is reached). This makes it possible to request a large number of nodes that continually pull in new jobs rather than submitting a large number of small jobs to the scheduler, which can be more efficient. In other words, don't be surprised if the Slurm job continues to run even when your submitted task has completed.

### Scaling Up

Now let's consider a more realistic scenario. Suppose we want to have a single Slurm job that reserves 8 nodes, and each `PythonApp` (e.g. VASP calculation) will run on 2 nodes (let's assume each node has 48 cores total, so that's a total of 96 cores for each calculation). Parsl will act as an orchestrator in the background of one of the nodes. Our config will now look like the following.

```python
n_parallel_calcs = 4 # Number of quacc calculations to run in parallel
n_nodes_per_calc = 2 # Number of nodes to reserve for each calculation
n_cores_per_node = 48 # Number of CPU cores per node
vasp_parallel_cmd = f"srun -N {n_nodes_per_calc} --ntasks={n_cores_per_node*n_nodes_per_calc} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores"

config = Config(
    max_idletime=300,
    executors=[
        HighThroughputExecutor(
            label="quacc_HTEX",
            max_workers=n_parallel_calcs,
            cores_per_worker=1e-6,
            provider=SlurmProvider(
                account="MyAccountName",
                nodes_per_block=n_nodes_per_calc*n_parallel_calcs,
                scheduler_options="#SBATCH -q debug -C cpu",
                worker_init=f"source activate quacc && module load vasp && export QUACC_VASP_PARALLEL_CMD={vasp_parallel_cmd}",
                walltime="00:10:00",
                launcher = SimpleLauncher(),
                cmd_timeout=120,
                init_blocks=0,
                min_blocks=1,
                max_blocks=1,
            ),
        )
    ],
)
```

In addition to some modified parameters, there are some new ones here too. We set `cores_per_worker` to a small value here so that the pilot job (e.g. the Parsl orchestrator) is allowed to be oversubscribed with scheduling processes. Setting `init_blocks`, `min_blocks`, and `max_blocks` like above ensures the right number of tasks are run.

```{seealso}
Dr. Logan Ward has a nice example on YouTube describing a very similar example [here](https://youtu.be/0V4Hs4kTyJs?t=398).
```

## Visualization

Parsl comes with a dashboard utility to visualize executed workflows when using the `HighThroughputExecutor`. Refer to the [Monitoring and Visualization](https://parsl.readthedocs.io/en/stable/userguide/monitoring.html#visualization) section of the Parsl documentation for details.

## Learn More

That ends the Parsl section of the documentation. If you want to learn more about Parsl, you can read the [Parsl Documentation](https://parsl.readthedocs.io/en/stable/#). Please refer to the [Parsl Slack Channel](http://parsl-project.org/support.html) for any Parsl-specific questions.
