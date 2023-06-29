# Using Quacc with FireWorks

## Introduction

[FireWorks](https://materialsproject.github.io/fireworks/) is a powerful software package to manage and execute complex workflows. It is best used in tandem with Jobflow because Jobflow comes with native support to convert a `Job` or `Flow` into a FireWorks `firework` or `workflow`, respectively.

```{hint}
Make sure you have completed the ["FireWorks Setup"](../../install/advanced/fireworks.md) instructions. Additionally, refer to the ["Using Quacc with Jobflow"](jobflow.md) section before reviewing the content below.
```

## Converting Between Jobflow and FireWorks

The {obj}`jobflow.managers.fireworks` module has all the tools you need to convert your Jobflow workflows to a format that is suitable for FireWorks.

### Converting a Job to a Firework

To convert a `Job` to a `firework` and add it to your launch pad:

```python
from fireworks import LaunchPad
from jobflow.managers.fireworks import job_to_firework

fw = job_to_firework(job)
lpad = LaunchPad.auto_load()
lpad.add_wf(fw)
```

### Converting a Flow to a Workflow

To convert a `Flow` to a `workflow` and add it to your launch pad:

```python
from fireworks import LaunchPad
from jobflow.managers.fireworks import flow_to_workflow

wf = flow_to_workflow(flow)
lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
```

## Dispatching Calculations

With a workflow added to your launch pad, on the desired machine of choice, you can run `qlaunch rapidfire --nlaunches <N>` (where `<N>` is the number of jobs to submit) in the command line to submit your workflows to the job scheduler. Running `qlaunch rapidfire -m <N>` will ensure that `<N>` jobs are always submitted to the queue.

## Monitoring the Launchpad

The easiest way to monitor the state of your launched FireWorks and workflows is through the GUI, which can be viewed with `lpad webgui`. To get the status of running fireworks from the command line, you can run `lpad get_fws -s RUNNING`. Other statuses can also be provided as well as individual FireWorks IDs.

To rerun a specific FireWork, one can use the `rerun_fws` command like so: `lpad rerun_fws -i <FWID>` where `<FWID>` is the FireWork ID. Similarly, one can rerun all fizzled jobs via `lpad rerun_fws -s FIZZLED`. More complicated Mongo-style queries can also be carried out. Cancelling a workflow can be done with `lpad delete_wflows -i <FWID>`.

Refer to the `lpad -h` help menu for more details.

## Continuous Job Submission

To ensure that jobs are continually submitted to the queue you can use `tmux` to preserve the job submission process even when the SSH session is terminated. For example, one can do the following to create a new `tmux` session called `launcher`:

```bash
module load tmux
tmux new -s launcher
```

To exit the `tmux` session, press `ctrl+b` followed by `d`. To re-enter the tmux session, run `tmux attach -t launcher`. Additional `tmux` commands can be found on the [tmux cheatsheet](https://tmuxcheatsheet.com/).

## Learn More

For details about how to assign different compute resources and submission options for each job in a workflow, refer to the [Setting where Jobs are Dispatched](https://materialsproject.github.io/jobflow/tutorials/8-fireworks.html#setting-where-jobs-are-dispatched) section in the Jobflow documentation.

For additional details on how to submit jobs to the queue that are in your launchpad, refer to the [FireWorks Documentation](https://materialsproject.github.io/fireworks/queue_tutorial.html#submit-a-job). Please refer to the [FireWorks MatSci Forum](https://matsci.org/c/fireworks/15) for FireWorks-specific questions.
