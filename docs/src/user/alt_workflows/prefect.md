# Using Quacc with Prefect

## Introduction

[Prefect](https://www.prefect.io/) is a workflow management system that is widely adopted in the data science industry. It can be used in place of Covalent, if preferred.

```{note}
For some minimal working examples of how to write your own Prefect workflows and how they compare to other workflow tools, refer to the [Worfklow Engine Comparison Guide](../comparison.md)
```

## Pre-Requisites

Make sure you completed the ["Prefect Setup"](../../install/alt_workflows/prefect.md) section of the documentation. Additionally, you should learn about the main Prefect concepts of a [`Flow`](https://docs.prefect.io/concepts/flows/) and a [`Task`](https://docs.prefect.io/concepts/tasks/), as described in the [Prefect Tutorial](https://docs.prefect.io/tutorial/).

## Examples

```{hint}
If you haven't logged into [Prefect Cloud](https://app.prefect.cloud/) yet, you should do so via `prefect cloud login`.
```

### Running a Simple Serial Workflow

We will now try running a simple workflow where we relax a bulk Cu structure using EMT and take the output of that calculation as the input to a follow-up static calculation with EMT.

```python
from prefect import flow, task
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job


# Define the workflow
@flow
def workflow(atoms):

    # Call Task 1
    future1 = task(relax_job).submit(atoms)

    # Call Task 2, which takes the output of Task 1 as input
    future2 = task(static_job).submit(future1.result()["atoms"])

    return future2.result()

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run the workflow with Prefect tracking
result = workflow(atoms)
print(result)
```

![Prefect UI](../../_static/user/prefect_tutorial1.jpg)

```{note}
We have used a short-hand notation here of `task(<function>)`. This is equivalent to using the `@task` decorator and definining a new function for each task. Calling `.submit()` enables concurrent execution of the tasks, which also requires the use of `.result()` to retrieve the output of the task.
```

![Prefect UI](../../_static/user/prefect_tutorial1.jpg)

### Running a Simple Parallel Workflow

Now let's consider a similar but nonetheless distinct example. Here, we will define a workflow where we will carry out two EMT structure relaxations, but the two jobs are not dependent on one another. In this example, Covalent will know that it can run the two jobs separately, and even if Job 1 were to fail, Job 2 would still progress.

```python
from prefect import flow
from ase.build import bulk, molecule
from quacc.recipes.emt.core import relax_job

# Define workflow
@flow
def workflow(atoms1, atoms2):

    # Define two independent relaxation jobs
    future1 = task(relax_job).submit(atoms1)
    future2 = task(relax_job).submit(atoms2)

    return {"result1": future1.result(), "result2": future2.result()}

# Define two Atoms objects
atoms1 = bulk("Cu")
atoms2 = molecule("N2")

# Run the workflow with Prefect tracking
result = workflow(atoms1, atoms2)
print(result)
```

![Prefect UI](../../_static/user/prefect_tutorial2.jpg)

### Running Workflows with Complex Connectivity

#### The Inefficient Way

For this example, let's consider a toy scenario where we wish to relax a bulk Cu structure, carve all possible slabs, and then run a new relaxation calculation on each slab. This is an example of a dynamic workflow.

In Quacc, there are two types of recipes: individual compute tasks with the suffix `_job` and pre-made multi-step workflows with the suffix `_flow`. Here, we are interested in importing a pre-made workflow. Refer to the example below:

```python
from prefect import task, flow
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow


@flow
def workflow(atoms):
    future1 = task(relax_job).submit(atoms)
    future2 = task(bulk_to_slabs_flow).submit(future1.result()["atoms"])

    return future2.result()

# Define the Atoms object
atoms = bulk("Cu")

# Run the workflow
result = workflow(atoms)
print(result)
```

When running a Covalent-based workflow like {obj}`.emt.slabs.bulk_to_slabs_flow` above, the entire function will run as a single compute task even though it is composed of several individual sub-tasks. If these sub-tasks are compute-intensive, this might not be the most efficient use of resources.

![Prefect UI](../../_static/user/prefect_tutorial3.jpg)

#### The Efficient Way

Quacc fully supports Prefect-based workflows to resolve this limitation. For example, the workflow above can be equivalently run as follows using the Prefect-specific {obj}`.emt.prefect.slabs.bulk_to_slabs_flow` workflow.

```python
from prefect import flow, task
from ase.build import bulk
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.prefect.slabs import bulk_to_slabs_flow

@flow
def workflow(atoms):
    future1 = task(relax_job).submit(atoms)
    result = bulk_to_slabs_flow(future1.result()["atoms"])

    return result

atoms = bulk("Cu")
result = workflow(atoms)
print(result)
```

In this example, all the individual tasks and sub-tasks are run as separate jobs, which is more efficient. By comparing {obj}`.emt.prefect.slabs.bulk_to_slabs_flow` with its Covalent counterpart {obj}`.emt.slabs.bulk_to_slabs_flow`, you can see that the two are extremely similar such that it is often straightforward to interconvert between the two. The `bulk_to_slabs_flow` used here is a Prefect `Flow` object, which is why we didn't need to wrap it with a `task()`. Since this is a `Flow` within a `Flow`, we call the inner flow a ["subflow"](https://docs.prefect.io/concepts/flows/?h=subflow#composing-flows).

![Prefect UI](../../_static/user/prefect_tutorial4.gif)

## Job Management

### Using a Prefect Agent

So far, we have dispatched calculations immediately upon calling them. However, in practice, it is often more useful to have a [Prefect agent](https://docs.prefect.io/concepts/work-pools/#agent-overview) running in the background that will continually poll for work to submit to the generated Dask cluster and job scheduler. This allows you to submit only a subset of workflows at a time, and the agent will automatically submit more jobs as the resources become available.

To run Prefect workflows with an agent, on the HPC environment where you wish to submit jobs, run `prefect agent start -p "quacc-pool"` to start a worker pool named "quacc-pool". Then submit your workflows as usual. It is best to run the agent on some perpetual resource like a login node or a dedicated workflow node.

### Defining a Task Runner

To modify where tasks are run, set the `task_runner` keyword argument of the corresponding `@flow` decorator. An example is shown below for setting up a `SLURMCluster` compatible with the NERSC Perlmutter machine. The jobs in this scenario would be submitted from a login node, and `prefect cloud login` should be run before submitting the workflow.

```python
from quacc.util.wflows import make_dask_cluster

n_jobs = 1 # Number of Slurm jobs
n_nodes = 1 # Number of nodes per Slurm job

cluster_params = {
    # Dask worker options
    "cores": 1,
    "memory": "4GB",
    "processes": 1,
    # SLURM options
    "shebang": "#!/bin/bash",
    "account": "MyAccountName",
    "walltime": "00:10:00",
    "job_mem": "0",
    "job_script_prologue": ["source ~/.bashrc", "conda activate quacc"],
    "job_directives_skip": ["-n", "--cpus-per-task"],
    "job_extra_directives": ["-q debug", f"-N {n_nodes}", "-C cpu"],
    "python": "python",
}

cluster = make_dask_cluster(cluster_params, n_jobs=n_jobs)
```

```{seealso}
Refer to the [Dask-Jobqueue Documentation](https://jobqueue.dask.org/en/latest/index.html) for the available keyword arguments to the Dask-generated clusters.
```

With this instantiated cluster object, you can set the task runner of a `Flow` as follows.

```python
from prefect_dask.task_runners import DaskTaskRunner

@flow(task_runner=DaskTaskRunner(cluster.scheduler_address))
def workflow(atoms):
    ...
```

Now, when the worklow is run from the login node, it will be submitted to the job scheduling system (Slurm in this case), and the results will be sent back to Prefect Cloud once completed. To modify an already imported `Flow` object, the `Flow.task_runner` attribute can also be modified directly.

## Learn More

That ends the Perfect section of the documentation. If you want to learn more about Perfect, you can read the [Prefect Documentation](https://docs.prefect.io/). Please refer to the [Prefect Slack Channel](https://www.prefect.io/slack/) and/or [Prefect Community Discourse](https://discourse.prefect.io/) page for any Prefect-specific questions.
