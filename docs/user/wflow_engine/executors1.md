# Deploying Calculations

In the previous examples, we have been running calculations on our local machine. However, in practice, you will probably want to run your calculations on one or more HPC machines. This section will describe how to set up your workflows to run on HPC machines using your desired workflow engine to scale up your calculations.

=== "Covalent"

    By default, Covalent will run all jobs on your local machine using the Dask backend. This is a parameter that you can control. For instance, Covalent offers many [executor plugins](https://docs.covalent.xyz/docs/plugin) that can be installed and used to interface with a wide range of HPC, cloud, and quantum devices.

    **Setting the Executor for the Flow**

    If you want to use the same executor for all the jobs in a workflow, you can pass the `executor` keyword argument to the `#!Python @flow` decorator.

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job


    @flow(executor="local")  # (1)!
    def workflow(atoms):
        result1 = relax_job(atoms)
        result2 = static_job(result1["atoms"])

        return result2


    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    1. This was merely for demonstration purposes. There is never really a need to use the "local" executor since the "dask" executor runs locally and is faster.

    **Setting Executors for Individual Jobs**

    The individual executor options for each job can be modified after they are imported as well.

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job


    @ct.electron(executor="local")
    def local_static_job(*args, **kwargs):
        return static_job(*args, **kwargs)


    @flow
    def workflow(atoms):
        output1 = relax_job(atoms)
        output2 = local_static_job(output1["atoms"])

        return output2


    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    **Configuring Executors**

    Refer to the [executor plugin documentation](https://docs.covalent.xyz/docs/plugin) for instructions on how to install and use the relevant plugins that allow Covalent to submit jobs on your desired machines. Most users of quacc will probably want to use the [`SlurmExecutor`](https://github.com/AgnostiqHQ/covalent-slurm-plugin), which is a plugin for Covalent that supports Slurm job scheduling system. If you are using a job scheduler environment but find the `covalent-slurm-plugin` does not suit your needs, you may wish to consider the [`covalent-hpc-plugin`](https://github.com/Quantum-Accelerators/covalent-hpc-plugin).

=== "Dask"

    A Dask cluster can be set up to be used with a queueing system like that found on most HPC machines. This is done via [Dask Jobqueue](https://jobqueue.dask.org/en/latest/index.html). Example configurations for various queuing systems can be found in the ["Example Deployments"](https://jobqueue.dask.org/en/latest/examples.html) section of the Dask Jobqueue documentation.

    !!! Note "Pilot Jobs"

        Unlike most other workflow engines, Dask is built for the [pilot job model](https://en.wikipedia.org/wiki/Pilot_job) where the allocated nodes continually pull in new tasks to run. This makes it possible to avoid submitting a large number of small jobs to the scheduler, which can be inefficient from a queuing perspective.

=== "Parsl"

    Out-of-the-box, Parsl will run on your local machine. However, in practice you will probably want to run your Parsl workflows on HPC machines.

    **Configuring Executors**

    To configure Parsl for the high-performance computing environment of your choice, refer to the [executor configuration page in the Parsl documentation](https://parsl.readthedocs.io/en/stable/userguide/configuring.html) for many examples. Additional details can be found in the ["Execution" section](https://parsl.readthedocs.io/en/stable/userguide/execution.html) of the Parsl documentation. Most users of quacc will probably want to the [`HighThroughputExecutor`](https://parsl.readthedocs.io/en/stable/stubs/parsl.executors.HighThroughputExecutor.html#parsl.executors.HighThroughputExecutor).

    !!! Note "Pilot Jobs"

        Unlike most other workflow engines, Parsl is built for the [pilot job model](https://en.wikipedia.org/wiki/Pilot_job) where the allocated nodes continually pull in new tasks to run. This makes it possible to avoid submitting a large number of small jobs to the scheduler, which can be inefficient from a queuing perspective.

=== "Prefect"

    To scale up calculations, read about the concept of a Prefect [task runner](https://docs.prefect.io/latest/concepts/task-runners/). By default, `quacc` automatically submits all `#!Python @job`-decorated functions to the specified task runner and so concurrency is achieved by default.

    To use Prefect in a job scheduler environment, you can create a [`DaskTaskRunner`](https://prefecthq.github.io/prefect-dask/usage_guide/) that can be used in conjunction with [dask-jobqueue](https://jobqueue.dask.org/en/latest). Example configurations for various queuing systems can be found in the ["Example Deployments"](https://jobqueue.dask.org/en/latest/examples.html) section of the `dask-jobqueue` documentation.

=== "Redun"

    Out-of-the-box, Redun will run on your local machine. However, in practice, you will probably want to specify a dedicated executor.

    To configure Redun for the high-performance computing environment of your choice, refer to the [executors](https://insitro.github.io/redun/executors.html) page in the Redun documentation.

=== "Jobflow"

    Out-of-the-box, Jobflow can be used to run on your local machine. You will, however, need a "manager" to run your workflows on HPC machines. The two recommended options are [jobflow-remote](https://github.com/Matgenix/jobflow-remote) and [FireWorks](https://github.com/materialsproject/fireworks), the latter of which is described here.

    **Setting Up Your `my_qadapter.yaml`**

    When you previously [set up Jobflow and FireWorks](../../install/wflow_engines.md), you created a `my_qadapter.yaml` file. It's now time to revisit that file and adjust the `pre_rocket` command with any modules or environment variables necessary for your calculations to run. Additionally, you will probably want to update the `nodes`, `walltime`, and related settings for your scheduler.

    **Converting Between Jobflow and FireWorks**

    The [`jobflow.managers.fireworks`](https://materialsproject.github.io/jobflow/jobflow.managers.html#module-jobflow.managers.fireworks) module has all the tools you need to convert your Jobflow workflows to a format that is suitable for FireWorks.

    **Converting a Job to a Firework**

    To convert a `Job` to a `firework` and add it to your launch pad:

    ```python
    from fireworks import LaunchPad
    from jobflow.managers.fireworks import job_to_firework

    fw = job_to_firework(job)
    lpad = LaunchPad.auto_load()
    lpad.add_wf(fw)
    ```

    **Converting a Flow to a Workflow**

    To convert a `Flow` to a `workflow` and add it to your launch pad:

    ```python
    from fireworks import LaunchPad
    from jobflow.managers.fireworks import flow_to_workflow

    wf = flow_to_workflow(flow)
    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf)
    ```

    **Setting Where Jobs are Dispatched**

    The `my_qadapter.yaml` file you made in the [installation instructions](../../install/install.md) specifies how FireWorks will submit jobs added to your launch pad. Additional details can be found in the [Jobflow Documentation](https://materialsproject.github.io/jobflow/tutorials/8-fireworks.html#setting-where-jobs-are-dispatched) for how to dynamically set where and how Jobflow `Job` and `Flow` objects can be dispatched.
