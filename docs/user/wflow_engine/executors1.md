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

    static_job.electron_object.executor = "local"

    @flow
    def workflow(atoms):
        output1 = relax_job(atoms)
        output2 = static_job(output1["atoms"])

        return output2


    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    ??? Tip "An Alternate Approach"

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

    !!! Tip "Example Configurations"

        Refer to the [executor plugin documentation](https://docs.covalent.xyz/docs/plugin) for instructions on how to install and use the relevant plugins that allow Covalent to submit jobs on your desired machines.

    Most users of quacc will probably want to use the [`HPCExecutor`](https://github.com/Quantum-Accelerators/covalent-hpc-plugin), which is a plugin for Covalent that supports Slurm, PBS, LSF, Flux, and more. For submitting jobs to a Slurm-based job scheduler from your local machine, an example `HPCExecutor` configuration might look like the following, which has been tested on Perlmutter at NERSC:

    ```python
    n_nodes = 2  # Number of nodes to reserve for each calculation
    n_cores_per_node = 48  # Number of CPU cores per node

    executor = ct.executor.HPCExecutor(
        # SSH credentials
        username="YourUserName",
        address="perlmutter-p1.nersc.gov",
        ssh_key_file="~/.ssh/nersc",
        cert_file="~/.ssh/nersc-cert.pub",  # (1)!
        # PSI/J parameters
        instance="slurm",
        resource_spec_kwargs={
            "node_count": n_nodes,
            "processes_per_node": n_cores_per_node,
        },  # (2)!
        job_attributes_kwargs={
            "duration": 10,  # minutes
            "project_name": "YourAccountName",
            "custom_attributes": {"slurm.constraint": "cpu", "slurm.qos": "debug"},
        },  # (3)!
        # Remote Python env parameters
        remote_conda_env="quacc",
        # Covalent parameters
        remote_workdir="$SCRATCH/quacc",  # (4)!
        create_unique_workdir=True,  # (5)!
        cleanup=False,  # (6)!
    )
    ```

    1. This a certificate file used to validate your SSH credentials. This is often not needed but is required at NERSC facilities due to the use of [`sshproxy`](https://docs.nersc.gov/connect/mfa/#sshproxy)-based multi-factor authentication.

    2. These are the resource specifications for the compute job, which are keyword arguments passed to PSI/J's [`ResourceSpecV1` class](https://exaworks.org/psij-python/docs/v/0.9.0/.generated/psij.html#psij.resource_spec.ResourceSpecV1).

    3. These are the job attributes that the job scheduler needs, which are keyword arguments passed to PSI/J's [`JobAttributes` class](https://exaworks.org/psij-python/docs/v/0.9.0/.generated/psij.html#psij.JobAttributes).

    4. This will tell Slurm where to `cd` and will specify where the calculations and results are stored. If you use this keyword argument, you should not explicitly specify the `RESULTS_DIR` quacc setting, which seeks to do largely the same thing.

    5. You generally want each quacc job to have the results stored in its own unique working directory to ensure files don't overwrite one another, so  `create_unique_workdir` should be set to `True`. This should not be used in combination with the `CREATE_UNIQUE_DIR` quacc setting, which seeks to do largely the same thing.

    6. For debugging purposes, it can be useful to keep all the temporary files. Once you're confident things work, you can omit the `cleanup` keyword argument.

    ??? Note "The SlurmExecutor"

        If you plan to use the dedicated [SlurmExecutor](https://docs.covalent.xyz/docs/user-documentation/api-reference/executors/slurm) developed by Covalent, an analogous example is included below:

        ```python
        n_nodes = 2
        n_cores_per_node = 48

        executor = ct.executor.SlurmExecutor(
            username="YourUserName",
            address="perlmutter-p1.nersc.gov",
            ssh_key_file="~/.ssh/nersc",
            cert_file="~/.ssh/nersc-cert.pub",
            conda_env="quacc",
            options={
                "nodes": f"{n_nodes}",
                "qos": "debug",
                "constraint": "cpu",
                "account": "YourAccountName",
                "job-name": "quacc",
                "time": "00:10:00",
            },
            remote_workdir="$SCRATCH/quacc",  # (1)!
            use_srun=False,  # (2)!
        )
        ```

        1. This will tell Slurm where to `cd` into and will specify where the calculations and results are stored. If you use this keyword argument, you should not explicitly specify the `RESULTS_DIR` quacc setting, which seeks to do largely the same thing.

        2.  The `SlurmExecutor` must have `use_srun=False` in order for ASE-based calculators to be launched appropriately.

=== "Dask"

    A Dask cluster can be set up to be used with a queueing system like that found on most HPC machines. This is done via [Dask Jobqueue](https://jobqueue.dask.org/en/latest/index.html). Example configurations for various queuing systems can be found in the ["Example Deployments"](https://jobqueue.dask.org/en/latest/examples.html) section of the documentation.

=== "Parsl"

    Out-of-the-box, Parsl will run on your local machine. However, in practice you will probably want to run your Parsl workflows on HPC machines.

    !!! Note "Pilot Jobs"

        Unlike most other workflow engines, Parsl is built for the [pilot job model](https://en.wikipedia.org/wiki/Pilot_job) where the allocated nodes continually pull in new tasks to run. This makes it possible to avoid submitting a large number of small jobs to the scheduler, which can be inefficient from a queuing perspective.

    **Configuring Executors**

    !!! Tip "Example Configurations"

        To configure Parsl for the high-performance computing environment of your choice, refer to the executor [Configuration page in the Parsl documentation](https://parsl.readthedocs.io/en/stable/userguide/configuring.html) for many examples.

    Let's imagine a scenario where we want to run a series of compute-intensive DFT calculations. Each DFT calculation requires 2 CPU nodes (each node having 48 cores). We want to run DFT calculations on 4 unique structures, all in parallel. We also want this to be done in a single Slurm allocation, meaning that this allocation must request 8 total nodes.

    An example [`HighThroughputExecutor`](https://parsl.readthedocs.io/en/stable/stubs/parsl.executors.HighThroughputExecutor.html#parsl.executors.HighThroughputExecutor) that will orchestrate jobs from the login node of NERSC's Perlmutter machine to carry out the above example is as follows:

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    max_slurm_jobs = 1  # Maximum number of Slurm jobs (blocks) to allocate
    n_calcs_per_job = 4  # Number of calculations to run in parallel (per block)
    n_nodes_per_calc = 2  # Number of nodes to reserve for each calculation

    config = Config(
        max_idletime=60,  # (1)!
        strategy="htex_auto_scale",  # (2)!
        executors=[
            HighThroughputExecutor(
                label="quacc_parsl",  # (3)!
                max_workers=n_calcs_per_job,  # (4)!
                cores_per_worker=1e-6,  # (5)!
                provider=SlurmProvider(
                    account="MyAccountName",
                    qos="debug",
                    constraint="cpu",
                    worker_init=f"source ~/.bashrc && conda activate quacc",  # (6)!
                    walltime="00:10:00",  # (7)!
                    nodes_per_block=n_nodes_per_calc * n_calcs_per_job,  # (8)!
                    init_blocks=0,  # (9)!
                    min_blocks=0,  # (10)!
                    max_blocks=max_slurm_jobs,  # (11)!
                    launcher=SimpleLauncher(),  # (12)!
                    cmd_timeout=120,  # (13)!
                ),
            )
        ],
    )

    parsl.load(config)
    ```

    1. The maximum amount of time (in seconds) to allow the executor to be idle before blocks (i.e. Slurm jobs) can potentially be shut down. Default is 120.

    2. Unique to the `HighThroughputExecutor`, this `strategy` will automatically scale the number of active blocks (i.e. Slurm jobs) up or down based on the number of tasks remaining. We set `max_blocks=1` here, so it can't scale up beyond 1 Slurm job, but it can scale down from 1 to 0 since `min_blocks=0`.

    3. This is just an arbitrary label for file I/O.

    4. Sets the maximum number of workers per block, which should generally be the number of tasks per block.

    5. This prevents the `HighThroughputExecutor` from reducing the number of workers if you request more workers than cores. It is [recommended](https://parsl.readthedocs.io/en/stable/userguide/mpi_apps.html#configuring-the-executor) for codes that run via MPI.

    6. Any additional `#SBATCH` options not captured elsewhere can be included here.

    7. Any commands to run before carrying out any of the Parsl tasks. This is useful for setting environment variables, activating a given Conda environment, and loading modules.

    7. The walltime for each block (i.e. Slurm job).

    8. The number of nodes that each block (i.e. Slurm job) should allocate.

    9. Sets the number of blocks (e.g. Slurm jobs) to provision during initialization of the workflow. We set this to a value of 0 so that there isn't a running Slurm job before any tasks have been submitted to Parsl.

    10. Sets the minimum number of blocks (e.g. Slurm jobs) to maintain during [elastic resource management](https://parsl.readthedocs.io/en/stable/userguide/execution.html#elasticity). We set this to 0 so that Slurm jobs aren't running when there are no remaining tasks.

    11. Sets the maximum number of active blocks (e.g. Slurm jobs) during [elastic resource management](https://parsl.readthedocs.io/en/stable/userguide/execution.html#elasticity). We set this to 1 here, but it can be increased to have multiple Slurm jobs running simultaneously. Raising `max_blocks` to a larger value will allow the "htex_auto_scale" strategy to upscale resources as needed.

    12. The type of Launcher to use. `SimpleLauncher()` must be used instead of the commonly used `SrunLauncher()` to allow quacc subprocesses to launch their own `srun` commands.

    13. The maximum time to wait (in seconds) for the job scheduler info to be retrieved/sent.

    **Practical Deployment**

    For debugging purposes or when running only a small numbers of jobs, it is simple enough to run the Parsl process from an interactive Jupyter Notebook or IPython kernel on the remote machine. However, for practical deployment and to ensure jobs are continually submitted to the queue even when the SSH session is terminated, you can run the Parsl orchestration process on a login node and maintain its state via a program like `tmux` or `screen`.

    For example, running `tmux new -s launcher` will create a new `tmux` session named `launcher`. To exit the `tmux` session while still preserving any running tasks on the login node, press `ctrl+b` followed by `d`. To re-enter the tmux session, run `tmux attach -t launcher`. Additional `tmux` commands can be found on the [tmux cheatsheet](https://tmuxcheatsheet.com/).

    **Multiple Executors**

    Parsl supports tying specific executors to a given `PythonApp`, as discussed in the [Multi-Executor section](https://parsl.readthedocs.io/en/stable/userguide/execution.html#multi-executor) of the Parsl documentation.

    ??? Note "Guide for NERSC Users"

        If you are a user of NERSC HPC resources, they have a [dedicated Parsl guide](https://docs.nersc.gov/jobs/workflow/parsl/) that is worth checking out.

=== "Prefect"

    To scale up calculations, read about the concept of a Prefect [task runner](https://docs.prefect.io/latest/concepts/task-runners/). By default, `quacc` automatically submits all `#!Python @job`-decorated functions to the specified task runner and so concurrency is achieved by default.

    To use Prefect in a job scheduler environment, you can create a [`DaskTaskRunner`](https://prefecthq.github.io/prefect-dask/usage_guide/) that can be used in conjunction with [dask-jobqueue](https://jobqueue.dask.org/en/latest). Example configurations for various queuing systems can be found in the ["Example Deployments"](https://jobqueue.dask.org/en/latest/examples.html) section of the `dask-jobqueue` documentation.

=== "Redun"

    Out-of-the-box, Redun will run on your local machine. However, in practice, you will probably want to specify a dedicated executor.

    !!! Tip "Example Configurations"

        To configure Redun for the high-performance computing environment of your choice, refer to the [executors](https://insitro.github.io/redun/executors.html) page in the Redun documentation.

=== "Jobflow"

    Out-of-the-box, Jobflow can be used to run on your local machine. You will, however, need a "manager" to run your workflows on HPC machines. The currently recommended manager for Jobflow is FireWorks, which is described here.

    **Setting Up Your `my_qadapter.yaml`**

    When you [set up Jobflow and FireWorks](../../install/wflow_engines.md), you created a `my_qadapter.yaml` file. It's now time to revisit that file and adjust the `pre_rocket` command with any modules or environment variables necessary for your calculations to run. Additionally, you will probably want to update the `nodes`, `walltime`, and related settings for your scheduler.

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

    **Dispatching Calculations**

    With a workflow added to your launch pad, on the desired machine of choice, you can run `qlaunch rapidfire --nlaunches <N>` (where `<N>` is the number of jobs to submit) in the command line to submit your workflows to the job scheduler. Running `qlaunch rapidfire -m <N>` will ensure that `<N>` jobs are always in the queue or running. To modify the order in which jobs are run, a priority can be set via `lpad set_priority <priority> -i <FWID>` where `<priority>` is a number.

    By default, `qlaunch` will launch compute jobs that each poll for a single FireWork to run. This means that more Slurm jobs may be submitted than there are jobs to run. To modify the behavior of `qlaunch` to only submit a Slurm job for each "READY" FireWork in the launchpad, use the `-r` ("reserved") flag.

    **Monitoring the Launchpad**

    The easiest way to monitor the state of your launched FireWorks and workflows is through the GUI, which can be viewed with `lpad webgui`. To get the status of running fireworks from the command line, you can run `lpad get_fws -s RUNNING`. Other statuses can also be provided as well as individual FireWorks IDs.

    To rerun a specific FireWork, one can use the `rerun_fws` command like so: `lpad rerun_fws -i <FWID>` where `<FWID>` is the FireWork ID. Similarly, one can rerun all fizzled jobs via `lpad rerun_fws -s FIZZLED`. More complicated Mongo-style queries can also be carried out. Cancelling a workflow can be done with `lpad delete_wflows -i <FWID>`.

    Refer to the `lpad -h` help menu for more details.

    **Continuous Job Submission**

    To ensure that jobs are continually submitted to the queue, you can use `tmux` to preserve the job submission process even when the SSH session is terminated. For example, running `tmux new -s launcher` will create a new `tmux` session named `launcher`. To exit the `tmux` session while still preserving any running tasks on the login node, press `ctrl+b` followed by `d`. To re-enter the tmux session, run `tmux attach -t launcher`. Additional `tmux` commands can be found on the [tmux cheatsheet](https://tmuxcheatsheet.com/).
