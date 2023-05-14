# Using QuAcc with Jobflow and FireWorks

## Jobflow

Make sure you completed the "Jobflow and Fireworks Setup" section of the installation instructions. Additionally, you should read the Jobflow [Quick Start](https://materialsproject.github.io/jobflow/tutorials/1-quickstart.html) to get a sense of how Jobflow works. Namely, you should understand the `Job` and `Flow` definitions, which describe individual compute tasks and workflows, respectively.

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

The key thing to note is that we need to transform the QuAcc recipe, which is a normal function, into a `jobflow` `Job` object. This can be done using the `@job` decorator and a new function definition or, more compactly, via `jf.job(<function>)`. We chose to run the job locally, but other workflow managers can be imported and used, as we discuss for FireWorks further below.

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

### Learn More

That ends the Jobflow section of the documentation. If you want to learn more about Jobflow, you can read the [Jobflow Documentation](https://materialsproject.github.io/jobflow/).

## FireWorks

[FireWorks](https://materialsproject.github.io/fireworks/) is a powerful software package to manage and execute complex workflows.

Jobflow comes with native support to convert a `Job` or `Flow` into a FireWorks `firework` or `workflow`, respectively.

To convert a `Job` to a `firework` and add it to your launch pad:

```python
from fireworks import LaunchPad
from jobflow.managers.fireworks import job_to_firework

fw = job_to_firework(job)
lpad = LaunchPad.auto_load()
lpad.add_wf(fw)
```

To convert a `Flow` to a `workflow` and add it to your launch pad:

```python
from fireworks import LaunchPad
from jobflow.managers.fireworks import flow_to_workflow

wf = flow_to_firework(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

With a workflow added to your launch pad, on the desired machine of choice, you can run `qlaunch rapidfire --nlaunches <N>` (where `<N>` is the number of jobs to submit) in the command line to submit your workflows to the job scheduler.

For additional FireWorks-related options in Jobflow, see the [Jobflow API](https://materialsproject.github.io/jobflow/jobflow.managers.html#module-jobflow.managers.fireworks). For documentation on how to submit jobs to the queue that are in your launchpad, refer to the [FireWorks Documentation](https://materialsproject.github.io/fireworks/queue_tutorial.html#submit-a-job).