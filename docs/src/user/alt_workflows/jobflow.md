# Using Quacc with Jobflow

## Introduction

[Jobflow](https://github.com/materialsproject/jobflow) is a program developed by the [Materials Project](https://materialsproject.org/) team to write computational workflows. It can be used in place of Covalent, if preferred.

```{note}
For some minimal working examples of how to write your own Jobflow workflows and how they compare to other tools, refer to the [Worfklow Engine Comparison Guide](alt_workflows/comparison.md).
```

## Pre-Requisites

Make sure you completed the ["Jobflow Setup"](../../install/alt_workflows/jobflow.md) section of the installation instructions. Additionally, you should read the Jobflow documentation's [Quick Start](https://materialsproject.github.io/jobflow/tutorials/1-quickstart.html) to get a sense of how Jobflow works. Namely, you should understand the `Job` and `Flow` definitions, which describe individual compute tasks and workflows, respectively.

```{seealso}
For a more detailed tutorial on how to use Jobflow, refer to the [Jobflow Tutorials](https://materialsproject.github.io/jobflow/tutorials) and [this helpful guide](https://github.com/JaGeo/Advanced_Jobflow_Tutorial) written by Dr. Janine George.
```

## Examples

### Running a Job

Here, we will try running a simple job where we carry out a static calculation on a bulk Cu structure using EMT. This example is shown below.

```python
import jobflow as jf
from ase.build import bulk
from quacc.recipes.emt.core import static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Define the compute job
job = jf.job(static_job)(atoms)

# Run the job locally
responses = jf.run_locally(job, create_folders=True)

# Get the result
result = responses[job.uuid][1].output
print(result)
```

The key thing to note is that we need to transform the Quacc recipe, which is a normal function, into a `Job` object. This can be done using the `@job` decorator and a new function definition or, more compactly, via `jf.job(<function>)`. We chose to run the job locally, but other workflow managers supported by Jobflow can be imported and used.

```{note}
Even though the Quacc recipes are defined as Covalent `Electron` objects via the `@ct.electron` decorator, this decorator will be ignored when using Jobflow.
```

### Running a Flow

Here, we will try running a simple workflow where we relax a bulk Cu structure using EMT and take the output of that calculation as the input to a follow-up static calculation with EMT. This example is shown below.

```python
import jobflow as jf
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Define the compute jobs
job1 = jf.job(relax_job)(atoms)
job2 = jf.job(static_job)(job1.output["atoms"])

# Define the workflow
workflow = jf.Flow([job1, job2])

# Run the workflow locally
responses = jf.run_locally(workflow, create_folders=True)

# Get the result
result = responses[job2.uuid][1].output
print(result)
```

Like before, we need to define the individual `Job` objects. Now though, we must stitch them together into a `Flow`, which can be easily achieved by passing them to the `jf.Flow()` constructor. The `Flow` object will automatically determine the order in which the jobs should be run based on the inputs and outputs of each job. In this case, it will know not to run `job2` until `job1` has completed.

### Running a Workflow with Complex Connectivity

#### The Inefficient Way

For this example, let's consider a toy scenario where we wish to relax a bulk Cu structure, carve all possible slabs, and then run a new relaxation calculation on each slab (with no static calculation at the end). This is an example of a dynamic workflow.

In Quacc, there are two types of recipes: individual compute tasks with the suffix `_job` and pre-made multi-step workflows with the suffix `_flow`. Here, we are interested in importing a pre-made workflow. Refer to the example below:

```python
import jobflow as jf
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow

# Define the Atoms object
atoms = bulk("Cu")

# Construct the Flow
job1 = jf.job(relax_job)(atoms)
job2 = jf.job(bulk_to_slabs_flow)(job1.output["atoms"])
workflow = jf.Flow([job1, job2])

# Run the workflow locally
responses = jf.run_locally(workflow, create_folders=True)

# Get the result
result = responses[job2.uuid][1].output
print(result)
```

We have imported the {obj}`.emt.slabs.bulk_to_slabs_flow` function, which takes an `Atoms` object along with several optional parameters. For demonstration purposes, we specify the `slab_static_electron=None` option to do a relaxation but disable the static calculation on each slab. All we have to do to define the workflow is stitch together the individual `@job` steps into a single `Flow` object.

#### The Efficient Way

Quacc fully supports Jobflow-based workflows to resolve this limitation. For example, the workflow above can be equivalently run as follows using the Jobflow-specific {obj}`.emt.jobflow.slabs.BulkToSlabsFlow` workflow:

```python
import jobflow as jf
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.jobflow.slabs import BulkToSlabsFlow

atoms = bulk("Cu")

@jf.job
def relax_func(atoms):
    return relax_job(atoms)

# Define the Atoms object
atoms = bulk("Cu")

# Construct the Flow
job1 = relax_func(atoms)
job2 = BulkToSlabsFlow().make(job1.output["atoms"])
workflow = jf.Flow([job1, job2])

# Run the workflow locally
jf.run_locally(workflow, store=STORE, create_folders=True, ensure_success=True)
```

In this example, all the individual tasks and sub-tasks are run as separate jobs, which is more efficient. By comparing {obj}`.emt.jobflow.slabs.BulkToSlabsFlow` with its Covalent counterpart {obj}`.emt.slabs.bulk_to_slabs_flow`, you can see that the two are extremely similar such that it is often straightforward to [interconvert](comparison.md) between the two. In the case of `BulkToSlabsFlow`, it actually returns a [`Response(replace)`](<https://materialsproject.github.io/jobflow/tutorials/5-dynamic-flows.html#The-Response(replace)-option>) object that dynamically replaces the `Flow` with several downstream jobs.

## Learn More

That ends the Jobflow section of the documentation. Continue to the ["Using Quacc with FireWorks"](fireworks.md) section to learn how to use Quacc with FireWorks, which is one of the manager options available with Jobflow to dispatch jobs in HPC environments.

If you want to learn more about Jobflow, you can read the [Jobflow Documentation](https://materialsproject.github.io/jobflow/). Please refer to the [Jobflow Discussions Board](https://github.com/materialsproject/jobflow/discussions) for Jobflow-specific questions.
