# Using Quacc with Jobflow

## Jobflow

Make sure you completed the ["Optional: Jobflow Setup"](../install/jobflow.md) and ["Optional: FireWorks Setup"](../install/fireworks.md) sections of the installation instructions. Additionally, you should read the Jobflow [Quick Start](https://materialsproject.github.io/jobflow/tutorials/1-quickstart.html) to get a sense of how Jobflow works. Namely, you should understand the `Job` and `Flow` definitions, which describe individual compute tasks and workflows, respectively.

### Example 1: Running a Job

Here, we will try running a simple job where we carry out a static calculation on a bulk Cu structure using EMT. This example is shown below.

```python
import jobflow as jf
from ase.build import bulk
from quacc.recipes.emt.core import static_job

atoms = bulk("Cu")

job = jf.job(static_job)(atoms)

responses = jf.run_locally(job, create_folders=True)

result = responses[job.uuid][1].output
print(result)
```

The key thing to note is that we need to transform the Quacc recipe, which is a normal function, into a `jobflow` `Job` object. This can be done using the `@job` decorator and a new function definition or, more compactly, via `jf.job(<function>)`. We chose to run the job locally, but other workflow managers can be imported and used, as we discuss for FireWorks further below.

### Example 2: Running a Flow

Here, we will try running a simple workflow where we relax a bulk Cu structure using EMT and take the output of that calculation as the input to a follow-up static calculation with EMT. This example is shown below.

```python
import jobflow as jf
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job

atoms = bulk("Cu")

job1 = jf.job(relax_job)(atoms)
job2 = jf.job(static_job)(job1.output["atoms"])

workflow = jf.Flow([job1, job2])

responses = jf.run_locally(workflow, create_folders=True)
result = responses[job2.uuid][1].output
print(result)
```

Like before, we need to define the individual `Job` objects. Now though, we must stitch them together into a `Flow`, which can be easily achieved by passing them to the `jf.Flow()` constructor. The `Flow` object will automatically determine the order in which the jobs should be run based on the inputs and outputs of each job. In this case, it will know not to run `job2` until `job1` has completed.

### Known Limitations

Jobflow cannot easily be used with Quacc recipes that involve classes (particularly those involving dynamic workflows) since they are structured for [Covalent](https://github.com/AgnostiqHQ/covalent). Nonetheless, QuAcc fully supports the development of Jobflow-specific workflows. Refer to the `quacc.recipes.emt.jobflow` module for an example.

### Learn More

That ends the Jobflow section of the documentation. If you want to learn more about Jobflow, you can read the [Jobflow Documentation](https://materialsproject.github.io/jobflow/). Continue to the ["Using Quacc with FireWorks"](fireworks.md) section to learn how to use Quacc with FireWorks.
