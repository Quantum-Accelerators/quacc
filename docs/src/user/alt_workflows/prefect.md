# Using Quacc with Prefect

## Introduction

## Examples

```{hint}
If you haven't logged into [Prefect Cloud](https://app.prefect.cloud/) yet, you may wish to do so via `prefect cloud login`.
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
    result1 = task(relax_job)(atoms)

    # Call Task 2, which takes the output of Task 1 as input
    result2 = task(static_job)(result1["atoms"])

    return result2

# Make an Atoms object of a bulk Cu structure
atoms = bulk("Cu")

# Run the workflow with Prefect tracking
result = workflow(atoms)
print(result)
```

![Prefect UI](../_static/user/prefect_tutorial.jpg)

```{note}
We have used a short-hand notation here of `task(<function>)`. This is equivalent to using the `@task` decorator and definining a new function for each task.
```

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
    result1 = task(relax_job)(atoms1)
    result2 = task(relax_job)(atoms2)

    return {"result1": result1, "result2": result2}

# Define two Atoms objects
atoms1 = bulk("Cu")
atoms2 = molecule("N2")

# Run the workflow with Prefect tracking
result = workflow(atoms1, atoms2)
print(result)
```

![Prefect UI](../_static/user/prefect_tutorial2.jpg)
