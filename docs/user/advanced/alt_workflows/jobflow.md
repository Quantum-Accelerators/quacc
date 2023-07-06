### Running a Job

Here, we will try running a simple job where we carry out a static calculation on a bulk Cu structure using EMT. This example is shown below.

```python
from jobflow import job, run_locally
from ase.build import bulk
from quacc.recipes.emt.core import static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Define the compute job
job1 = job(static_job)(atoms)

# Run the job locally
responses = run_locally(job1, create_folders=True)

# Get the result
result = responses[job1.uuid][1].output
print(result)
```

The key thing to note is that we need to transform the quacc recipe, which is a normal function, into a `Job` object. This can be done using the `@job` decorator and a new function definition or, more compactly, via `job(<function>)`. We chose to run the job locally, but other workflow managers supported by Jobflow can be imported and used.

```{note}
Even though the quacc recipes are defined as Covalent `Electron` objects via the `@ct.electron` decorator, this decorator will be ignored when using Jobflow.
```

### Running a Flow

Here, we will try running a simple workflow where we relax a bulk Cu structure using EMT and take the output of that calculation as the input to a follow-up static calculation with EMT. This example is shown below.

```python
from jobflow import job, Flow, run_locally
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Define the compute jobs
job1 = job(relax_job)(atoms)
job2 = job(static_job)(job1.output)

# Define the workflow
workflow = Flow([job1, job2])

# Run the workflow locally
responses = run_locally(workflow, create_folders=True)

# Get the result
result = responses[job2.uuid][1].output
print(result)
```

Like before, we need to define the individual `Job` objects. Now though, we must stitch them together into a `Flow`, which can be easily achieved by passing them to the `Flow()` constructor. The `Flow` object will automatically determine the order in which the jobs should be run based on the inputs and outputs of each job. In this case, it will know not to run `job2` until `job1` has completed.

### Running a Workflow with Complex Connectivity

#### The Inefficient Way

For this example, let's consider a toy scenario where we wish to relax a bulk Cu structure, carve all possible slabs, and then run a new relaxation calculation on each slab (with no static calculation at the end). This is an example of a dynamic workflow.

In quacc, there are two types of recipes: individual compute tasks with the suffix `_job` and pre-made multi-step workflows with the suffix `_flow`. Here, we are interested in importing a pre-made workflow. Refer to the example below:

```python
from jobflow immport job, Flow, run_locally
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow

# Define the Atoms object
atoms = bulk("Cu")

# Construct the Flow
job1 = job(relax_job)(atoms)
job2 = job(bulk_to_slabs_flow)(job1.output)
workflow = Flow([job1, job2])

# Run the workflow locally
responses = run_locally(workflow, create_folders=True)

# Get the result
result = responses[job2.uuid][1].output
print(result)
```

We have imported the `.emt.slabs.bulk_to_slabs_flow` function, which takes an `Atoms` object along with several optional parameters. For demonstration purposes, we specify the `slab_static_electron=None` option to do a relaxation but disable the static calculation on each slab. All we have to do to define the workflow is stitch together the individual `@job` steps into a single `Flow` object.

#### The Efficient Way

Quacc fully supports Jobflow-based workflows to resolve this limitation. For example, the workflow above can be equivalently run as follows using the Jobflow-specific `.emt.jobflow.slabs.bulk_to_slabs_flow` workflow:

```python
from jobflow import job, Flow, run_locally
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.jobflow.slabs import bulk_to_slabs_flow

# Define the Atoms object
atoms = bulk("Cu")

# Construct the Flow
job1 = job(relax_job)(atoms)
job2 = job(bulk_to_slabs_flow)(job1.output)
workflow = Flow([job1, job2])

# Run the workflow locally
run_locally(workflow, create_folders=True)
```

In this example, all the individual tasks and sub-tasks are run as separate jobs, which is more efficient. By comparing `.emt.jobflow.slabs.bulk_to_slabs_flow` with its Covalent counterpart `.emt.slabs.bulk_to_slabs_flow`, you can see that the two are extremely similar such that it is often straightforward to [interconvert](comparison.md) between the two. In the case of `bulk_to_slabs_flow`, it actually returns a [`Response(replace)`](<https://materialsproject.github.io/jobflow/tutorials/5-dynamic-flows.html#The-Response(replace)-option>) object that dynamically replaces the `Flow` with several downstream jobs.
