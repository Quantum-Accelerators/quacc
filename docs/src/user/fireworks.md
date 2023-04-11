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
