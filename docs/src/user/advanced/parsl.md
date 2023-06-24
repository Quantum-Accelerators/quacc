# Using Quacc with Parsl

## Introduction

[Parsl](https://github.com/Parsl/parsl) is a Python program developed at Argonne National Laboratory and the University of Chicago to easily write parallel workflows that can be dispatched on distributed compute resources. Like Jobflow+FireWorks, it can be used in place of Covalent, if preferred.

Make sure you completed the ["Parsl Setup"](../../install/advanced/parsl.md) section of the installation instructions. Additionally, you should read the Parsl documentation's ["Quick Start"](https://parsl.readthedocs.io/en/stable/quickstart.html) to get a sense of how Parsl works. Namely, you should understand the concept of a `@python_app` and `@join_app`, which describe individual compute tasks and dynamic job tasks, respectively.

```{seealso}
For a more detailed tutorial on how to use Parsl, refer to the ["Parsl Tutorial"](https://parsl.readthedocs.io/en/stable/1-parsl-introduction.html) and the even more detailed ["Parsl User Guide"](https://parsl.readthedocs.io/en/stable/userguide/index.html).
```

## Examples

### Running a Simple Serial Workflow

```{hint}
If you haven't loaded your Parsl config, you must do that first so Parsl can construct the job dependency graph. For testing purposes, you simply can run `import parsl` followed by `parsl.load()` before starting the examples below, which will enable jobs to run on your local machine.
```

We will first try running a simple workflow where we relax a bulk Cu structure using EMT and take the output of that calculation as the input to a follow-up static calculation with EMT.

```python
import parsl
from parsl import python_app
from ase.build import bulk

# Define the Python apps
@python_app
def relax_app(atoms):

    # All dependencies must be inside the Python app
    from quacc.recipes.emt.core import relax_job

    return relax_job(atoms)

@python_app
def static_app(atoms):

    # All dependencies must be inside the Python app
    from quacc.recipes.emt.core import static_job

    return static_job(atoms)

# Define the workflow
def workflow(atoms):

    # Call Job 1
    future1 = relax_app(atoms)

    # Call Job 2, which takes the output of Job 1 as input
    future2 = static_app(future1.result()["atoms"])

    return future2

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run the workflow
wf_future = workflow(atoms)
print(wf_future.result())
```

You can see that it is quite trivial to set up a Parsl workflow using the recipes within Quacc. We define the full workflow as a function that stitches together the individual `@python_app` workflow steps.

```{note}
The use of `.result()` serves to block any further calculations from running until it is resolved. Calling `.result()` also returns the function output as opposed to the `AppFuture` object.
```

```{warning}
Don't call `.result()` in a `return` statement. It will not block like you might naively expect it to.
```

### Running a Simple Parallel Workflow

Now let's consider a similar but nonetheless distinct example. Here, we will define a workflow where we will carry out two EMT structure relaxations, but the two jobs are not dependent on one another. In this example, Parsl will know that it can run the two jobs in parallel, and even if Job 1 were to fail, Job 2 would still progress.

```python
from parsl import python_app
from ase.build import bulk, molecule

# Define the Python app
@python_app
def relax_app(atoms):

    # All dependencies must be inside the Python app
    from quacc.recipes.emt.core import relax_job

    return relax_job(atoms)

# Define workflow
def workflow(atoms1, atoms2):

    # Define two independent relaxation jobs
    future1 = relax_app(atoms1)
    future2 = relax_app(atoms2)

    return future1, future2

# Define two Atoms objects
atoms1 = bulk("Cu")
atoms2 = molecule("N2")

# Run the workflow
future1, future2 = workflow(atoms1, atoms2)
print(future1.result(), future2.result())
```

If you monitor the output, you'll notice that the two jobs are being run in parallel, whereas they would be run sequentially if you were not using Parsl.

### Running Workflows with Complex Connectivity

For this example, let's consider a toy scenario where we wish to relax a bulk Cu structure, carve all possible slabs, and then run a new relaxation calculation on each slab (with no static calculation at the end).

In Quacc, there are two types of recipes: individual compute tasks with the suffix `_job` and pre-made multi-step workflows with the suffix `_flow`. Here, we are interested in importing a pre-made workflow. Refer to the example below:

```python
from parsl import python_app
from ase.build import bulk

@python_app
def relax_app(atoms):

    # All dependencies must be inside the Python app
    from quacc.recipes.emt.core import relax_job

    return relax_job(atoms)

@python_app
def bulk_to_slabs_app(atoms):

    # All dependencies must be inside the Python app
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

    return bulk_to_slabs_flow(atoms, slab_static_electron=None)

def workflow(atoms):
    future1 = relax_app(atoms)
    future2 = bulk_to_slabs_app(future1.result()["atoms"])

    return future2

# Define the Atoms object
atoms = bulk("Cu")

# Run the workflow
wf_future = workflow(atoms)
print(wf_future.result())
```

When running a Covalent-based workflow like {obj}`.emt.slabs.bulk_to_slabs_flow` above, the entire function will run as a single compute task even though it is composed of several individual sub-tasks. If these sub-tasks are compute-intensive, this might not be the most efficient use of resources.

Quacc fully supports the development of Parsl-based workflows to resolve this limitation. For example, the workflow above can be equivalently run as follows using the Parsl-specific {obj}`.emt.parsl.slabs.bulk_to_slabs_flow` workflow:

```python
from parsl import python_app
from ase.build import bulk
from quacc.recipes.emt.parsl.slabs import bulk_to_slabs_flow

@python_app
def relax_app(atoms):

    from quacc.recipes.emt.core import relax_job

    return relax_job(atoms)

atoms = bulk("Cu")

relax_future = relax_app(atoms)

wf_future = bulk_to_slabs_flow(relax_future.result()["atoms"], slab_static_app=None)
print(wf_future.result())
```

In this example, all the individual tasks and sub-tasks are run as separate jobs, which is more efficient. By comparing {obj}`.emt.parsl.slabs.bulk_to_slabs_flow` with its Covalent counterpart {obj}`.emt.slabs.bulk_to_slabs_flow`, you can see that the two are extremely similar such that it is often straightforward to interconvert between the two.

```{note}
We didn't need to wrap `bulk_to_slabs_flow` with a decorator because it is defined in Quacc as a `@join_app` (similar to a `@python_app` for dynamic workflow steps) that itself returns an `AppFuture`. This is also why we call `.result()` on it.
```

## Visualization

Parsl comes with a web dashboard utility to visualize executed workflows. Refer to the [Monitoring and Visualization](https://parsl.readthedocs.io/en/stable/userguide/monitoring.html#visualization) section of the Parsl documentation for details.

## Setting Executors

Out-of-the-box, Parsl will run on your local machine. However, in practice you will probably want to run your Parsl workflows on HPC machines.

```{note}
If you are just starting out, try running some test calculations locally first. Then come back and set up the relevant configuration files for your desired machines.
```

### Configuring Executors

To configure Parsl for the high-performance computing environment of your choice, refer to the executor [Configuration](https://parsl.readthedocs.io/en/stable/userguide/configuring.html) page in the Parsl documentation.

For [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/), example `HighThroughputExecutor` configurations can be found in the [NERSC Documentation](https://docs.nersc.gov/jobs/workflow/parsl/). A simple one is reproduced below that allows for job submission from the login node:

```python
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SimpleLauncher
from parsl.providers import SlurmProvider

config = Config(
    executors=[
        HighThroughputExecutor(
            label="quacc",
            max_workers=1,
            provider=SlurmProvider(
                partition="debug",
                account="MyAccountName",
                nodes_per_block=1,
                scheduler_options="#SBATCH -C cpu",
                worker_init="source activate quacc",
                walltime="00:10:00",
                launcher = SimpleLauncher()
            ),
        )
    ]
)
```

```{warning}
This example still needs to be tested and will likely need to be updated slightly.
```

### Swapping Executor Configurations

```{warning}
This section on how to dynamically set executor options is still a work in progress.
```

## Learn More

That ends the Parsl section of the documentation. If you want to learn more about Parsl, you can read the [Parsl Documentation](https://parsl.readthedocs.io/en/stable/#). Please refer to the [Parsl Slack Channel](http://parsl-project.org/support.html) for any Parsl-specific questions.
