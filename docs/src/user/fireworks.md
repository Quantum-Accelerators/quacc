# Using QuAcc with Jobflow and FireWorks

## Jobflow

Make sure you completed the "Jobflow and Fireworks Setup" section of the installation instructions. Additionally, you should read the Jobflow [Quick Start](https://materialsproject.github.io/jobflow/tutorials/1-quickstart.html) to get a sense of how Jobflow works. Namely, you should understand the `Job` and `Flow` definitions, which describe individual compute tasks and workflows, respectively.

### Example 1: Running a Job

```python
import jobflow as jf
from ase.build import bulk
from jobflow.managers.local import run_locally
from quacc.recipes.emt.core import static_job

atoms = bulk("Cu")

job = jf.job(static_job)(atoms)

responses = run_locally(job)

result = responses[job.uuid][1].output
print(result)
```

### Example 2: Running a Flow

```python
import jobflow as jf
from ase.build import bulk
from jobflow.managers.local import run_locally
from quacc.recipes.emt.core import relax_job, static_job

atoms = bulk("Cu")

job1 = jf.job(relax_job)(atoms)
job2 = jf.job(static_job)(job1.output["atoms"])

workflow = jf.Flow([job1, job2])

responses = run_locally(workflow)
result = responses[job2.uuid][1].output
print(result)
```

## FireWorks

[FireWorks](https://materialsproject.github.io/fireworks/) is a powerful software package to manage and execute complex workflows.

Jobflow comes with native support to convert a `Job` or `Flow` into a FireWorks `firework` or `workflow`, respectively.

To convert a `Job` to a `firework` and add it to your launch pad:

```python
from fireworks import LaunchPad
from jobflow.managers.fireworks import job_to_firework

fw = job_to_firework(job)
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

To convert a `Flow` to a `workflow` and add it to your launch pad:

```python
from fireworks import LaunchPad
from jobflow.managers.fireworks import flow_to_workflow

wf = flow_to_firework(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

For additional options, see the [Jobflow API](https://materialsproject.github.io/jobflow/jobflow.managers.html#module-jobflow.managers.fireworks). For documentation on how to submit jobs to the queue that are in your launchpad, refer to the [FireWorks documentation](https://materialsproject.github.io/fireworks/queue_tutorial.html#submit-a-job).