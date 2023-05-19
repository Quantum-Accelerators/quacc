# Using Quacc with FireWorks

## Introduction

[FireWorks](https://materialsproject.github.io/fireworks/) is a powerful software package to manage and execute complex workflows. It is best used in tandem with Jobflow because Jobflow comes with native support to convert a `Job` or `Flow` into a FireWorks `firework` or `workflow`, respectively.

```{hint}
Make sure you have completed the ["FireWorks Setup"](../../install/advanced/fireworks.md) instructions. Additionally, refer to the ["Using Quacc with Jobflow"](jobflow.md) section before reviewing the content below.
```

## Converting a Job to a Firework

To convert a `Job` to a `firework` and add it to your launch pad:

```python
from fireworks import LaunchPad
from jobflow.managers.fireworks import job_to_firework

fw = job_to_firework(job)
lpad = LaunchPad.auto_load()
lpad.add_wf(fw)
```

## Converting a Flow to a Workflow

To convert a `Flow` to a `workflow` and add it to your launch pad:

```python
from fireworks import LaunchPad
from jobflow.managers.fireworks import flow_to_workflow

wf = flow_to_firework(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

## Dispatching Calculations

With a workflow added to your launch pad, on the desired machine of choice, you can run `qlaunch rapidfire --nlaunches <N>` (where `<N>` is the number of jobs to submit) in the command line to submit your workflows to the job scheduler.

## Read More

For additional FireWorks-related options in Jobflow, see the [Jobflow API](https://materialsproject.github.io/jobflow/jobflow.managers.html#module-jobflow.managers.fireworks). For documentation on how to submit jobs to the queue that are in your launchpad, refer to the [FireWorks Documentation](https://materialsproject.github.io/fireworks/queue_tutorial.html#submit-a-job).
