# Using Quacc with Jobflow

## Introduction

[Jobflow](https://github.com/materialsproject/jobflow) is a program developed by the [Materials Project](https://materialsproject.org/) team to write computational workflows. It can be used in place of Covalent, if preferred.

Make sure you completed the ["Jobflow Setup"](../../install/advanced/jobflow.md) section of the installation instructions. Additionally, you should read the Jobflow documentation's [Quick Start](https://materialsproject.github.io/jobflow/tutorials/1-quickstart.html) to get a sense of how Jobflow works. Namely, you should understand the `Job` and `Flow` definitions, which describe individual compute tasks and workflows, respectively.

```{note}
For a more detailed tutorial on how to use Jobflow, refer to the [Jobflow Tutorials](https://materialsproject.github.io/jobflow/tutorials) and [this helpful guide](https://github.com/JaGeo/Advanced_Jobflow_Tutorial) written by Dr. Janine George.
```

### Example 1: Running a Job

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

### Example 2: Running a Flow

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

### Known Limitations

Jobflow cannot be used with Quacc recipes that involve classes since they are, by default, structured for [Covalent](https://github.com/AgnostiqHQ/covalent). To address this, we fully support the development of Jobflow-specific workflows that are mirrors of their Covalent counterparts. Refer to the {obj}`quacc.recipes.emt.jobflow.slabs` module for a representative example that can be compared against the Covalent version at {obj}`quacc.recipes.emt.slabs`.

### Learn More

That ends the Jobflow section of the documentation. If you want to learn more about Jobflow, you can read the [Jobflow Documentation](https://materialsproject.github.io/jobflow/). Continue to the ["Using Quacc with FireWorks"](fireworks.md) section to learn how to use Quacc with FireWorks.
