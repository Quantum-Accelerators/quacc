# Using Quacc with Parsl

## Introduction

[Parsl](https://github.com/Parsl/parsl) is a Python program developed at Argonne National Laboratory and the University of Chicago to easily write parallel workflows that can be dispatched on distributed compute resources. Like Jobflow+FireWorks, it can be used in place of Covalent, if preferred.

Make sure you completed the ["Parsl Setup"](../../install/advanced/parsl.md) section of the installation instructions. Additionally, you should read the Parsl documentation's [Quick Start](https://parsl.readthedocs.io/en/stable/quickstart.html) to get a sense of how Parsl works. Namely, you should understand the concept of a `@python_app` and `@join_app`, which describe individual compute tasks and dynamic job tasks, respectively.

```{note}
For a more detailed tutorial on how to use Parsl, refer to the [Parsl Tutorial](https://parsl.readthedocs.io/en/stable/1-parsl-introduction.html) and the even more detailed [Parsl user guide](https://parsl.readthedocs.io/en/stable/userguide/index.html).
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

```{note}
Note that the use of `.result()` serves to block any further calculations from running until it is resolved. Calling `.result()` also returns the function output as opposed to the `AppFuture` object.
```

```{hint}
Don't call `.result()` in a `return` statement. It will not block like you might naively expect it to.
```

You can see that it is quite trivial to set up a Parsl workflow using the recipes within Quacc. We define the full workflow as a simple function that stitches together the individual `@python_app` workflow steps.

### Running a Simple Parallel Workflow

Now let's consider a similar but nonetheless distinct example. Here, we will define a workflow where we will carry out two EMT structure relaxations, but the two jobs are not dependent on one another. In this example, Parsl will know that it can run the two jobs in parallel, and even if Job 1 were to fail, Job 2 would still progress.

```python
import parsl
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
print(wf_result["result1"].result(), wf_result["result2"].result())
```

### Running Workflows with Complex Connectivity

For this example, let's consider a toy scenario where we wish to relax a bulk Cu structure, carve all possible slabs, and then run a new relaxation calculation on each slab.

In Quacc, there are two types of recipes: 1) individual compute tasks with the suffix `_job`; pre-made multi-step workflows with the suffix `_flow`. Here, we are interested in importing a pre-made workflow. Refer to the example below:

```python
from ase.build import bulk
from quacc.recipes.emt.parsl.slabs import bulk_to_slabs_flow

wf_result = bulk_to_slabs_flow(bulk("Cu"), slab_static_app=None)
print(wf_result)
```

```{note}
We have called `.result()` here because `bulk_to_slabs_flow` is a `@join_app` (similar to a `@python_app` for dynamic workflow steps) that returns an `AppFuture`.
```

We have imported the {obj}`.emt.slabs.parsl.bulk_to_slabs_flow` class, which is supplied an `Atoms` object. Here, for demonstration purposes, we specify the `slab_static_app=None` option to do a relaxation but disable the static calculation on each slab.

```{hint}
If you are interested in rewriting a Covalent workflow into Parsl, it is often relatively straightforward. Compare {obj}`quacc.recipes.emt.slabs` and {obj}`quacc.recipes.emt.slabs.parsl` for the key differences.
```

## Visualization

Parsl comes with a web dashboard utility to visualize executed workflows. Refer to the [Monitoring and Visualization](https://parsl.readthedocs.io/en/stable/userguide/monitoring.html#visualization) section of the Parsl documentation for details.

## Setting Executors

```{note}
If you are just starting out, try running some test calculations locally first. Then come back and set up the relevant configuration files for your desired machines.
```

Out-of-the-box, Parsl will run on your local machine. TODO!!!

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

## Learn More

That ends the Parsl section of the documentation. If you want to learn more about Parsl, you can read the [Parsl Documentation](https://parsl.readthedocs.io/en/stable/#). Please refer to the [Parsl Slack Channel](http://parsl-project.org/support.html) for any Parsl-specific questions.
```
