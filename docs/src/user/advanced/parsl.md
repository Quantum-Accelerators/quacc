# Using Quacc with Parsl

## Introduction

[Parsl](https://github.com/Parsl/parsl) is a Python program developed at Argonne National Laboratory and the University of Chicago to easily write parallel workflows that can be dispatched on distributed compute resources. Like Jobflow+FireWorks, it can be used in place of Covalent, if preferred.

Make sure you completed the ["Parsl Setup"](../../install/advanced/parsl.md) section of the installation instructions. Additionally, you should read the Parsl documentation's [Quick Start](https://parsl.readthedocs.io/en/stable/quickstart.html) to get a sense of how Parsl works. Namely, you should understand the concept of a `@python_app`, which describes individual compute tasks.

```{note}
For a more detailed tutorial on how to use Parsl, refer to the [Parsl Tutorial](https://parsl.readthedocs.io/en/stable/1-parsl-introduction.html) and the even more detailed [Parsl user guide](https://parsl.readthedocs.io/en/stable/userguide/index.html).
```

## Running a Simple Serial Workflow

```{hint}
If you haven't loaded your Parsl config, you must do that first so Parsl can construct the job dependency graph. For testing purposes, you simply can run `import parsl` followed by `parsl.load()` before starting the examples below`.
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

# Define the Workflow
def workflow(atoms):

    # Call Job 1
    future1 = relax_app(atoms)

    # Call Job 2, which takes the output of Job 1 as input
    # Note the use of .result(), which blocks until the result is ready
    future2 = static_app(future1.result()["atoms"])

    return future2.result()

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run the workflow
wf_result = workflow(atoms)
print(wf_result)
```

You can see that it is quite trivial to set up a workflow using the recipes within Quacc. We define the full workflow as a simple function that stitches together the individual `@python_app` workflow steps.

## Running a Simple Parallel Workflow

Now let's consider a similar but nonetheless distinct example. Here, we will define a workflow where we will carry out two EMT structure relaxations, but the two jobs are not dependent on one another. In this example, Parsl will know that it can run the two jobs in parallel, and even if Job 1 were to fail, Job 2 would still progress.

```python
import parsl
from parsl import python_app
from ase.build import bulk, molecule

# Define the Python apps
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

    return {"result1": future1.result(), "result2": future2.result()}

# Define two Atoms objects
atoms1 = bulk("Cu")
atoms2 = molecule("N2")

# Run the workflow
wf_result = workflow(atoms1, atoms2)
print(wf_result)
```

## Running Workflows with Complex Connectivity

For this example, let's consider a toy scenario where we wish to relax a bulk Cu structure, carve all possible slabs, and then run a new relaxation calculation on each slab.

In Quacc, there are two types of recipes: 1) individual compute tasks that are functions; 2) workflows that are classes. Here, we are interested in importing a workflow, so it will be instantiated slightly differently from the prior examples. See the example below:

```python
import parsl
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
    from quacc.recipes.emt.slabs import BulkToSlabsFlow

    return BulkToSlabsFlow(slab_static_electron=None).run(atoms)

def workflow(atoms):
    future1 = relax_app(atoms)
    future2 = bulk_to_slabs_app(future1.result()["atoms"])

    return future2.result()

# Define the Atoms object
atoms = bulk("Cu")

# Run the workflow
wf_result = workflow(atoms)
print(wf_result)
```

We have imported the {obj}`.emt.slabs.BulkToSlabsFlow` class, which is instantiated with optional parameters and is applied to an `Atoms` object. Here, for demonstration purposes, we specify the `slab_static_electron=None` option to do a relaxation but disable the static calculation on each slab. All we have to do to define the workflow is stitch together the individual `@python_app` steps into a single function.

### Known Limitations

When running a Covalent-based class like {obj}`.emt.slabs.BulkToSlabsFlow` in the previous example, the entire class will run as a single compute task even though it is composed of several individual sub-tasks. If these sub-tasks are compute-intensive, this might not be the most efficient use of resources.

To address this, you can draw inspiration from the Covalent-based classes to design your own workflows tailored to Parsl. After all, it only requires you to stitch together the individual `@python_app` steps into a single function.

```{seealso}
For details on how to write your own dynamic workflows in Parsl, refer to the `@join_app` section of the [Parsl documentation](https://parsl.readthedocs.io/en/stable/1-parsl-introduction.html#Examples).
```

If you wish to construct Parsl-specific workflows that are mirrors of their Covalent counterparts, this is fully supported by Quacc.

### Visualization

Parsl comes with a web dashboard utility to visualize executed workflows. Refer to the [Monitoring and Visualization](https://parsl.readthedocs.io/en/stable/userguide/monitoring.html#visualization) section of the Parsl documentation for details.
