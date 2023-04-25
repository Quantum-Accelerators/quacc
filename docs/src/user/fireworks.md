# Using QuAcc with FireWorks

[FireWorks](https://materialsproject.github.io/fireworks/) is a powerful software package to manage and execute complex workflows. It is particularly useful for high-throughput computational screening efforts to efficiently manage large numbers of compute jobs.

Jobflow comes with native support to convert a `job` or `flow` into a FireWorks `firework` or `workflow`, respectively.

To convert a `job` to a `firework` and add it to your launch pad:

```python
from fireworks import LaunchPad
from jobflow.managers.fireworks import job_to_firework

fw = job_to_firework(job)
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

To convert a `flow` to a `workflow` and add it to your launch pad:

```python
from fireworks import LaunchPad
from jobflow.managers.fireworks import flow_to_workflow

wf = flow_to_firework(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

For additional options, see the [Jobflow API](https://materialsproject.github.io/jobflow/jobflow.managers.html#module-jobflow.managers.fireworks).

## Jobflow+Fireworks

Show fireworks specifically.

### Example 1: Running a Single Function

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

### Example 2: Running a Simple Serial Workflow

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

### Example 3: Running a Simple Parallel Workflow

```python
import jobflow as jf
from ase.build import bulk, molecule
from quacc.recipes.emt.core import static_job
from jobflow.managers.local import run_locally

atoms1 = bulk("Cu")
atoms2 = molecule("N2")
job1 = jf.job(static_job)(atoms1)
job2 = jf.job(static_job)(atoms2)
responses1 = run_locally(job1)
responses2 = run_locally(job2)

print([result1, result2])
```
