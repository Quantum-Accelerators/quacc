# Deploying Calculations

In the previous examples, we have been running calculations on our local machine. However, in practice, you will probably want to run your calculations on one or more HPC machines. This section will describe how to set up your workflows to run on HPC machines using your desired workflow engine to scale up your calculations.

=== "Covalent"

    By default, Covalent will run all `Electron` tasks on your local machine using the [`DaskExecutor`](https://docs.covalent.xyz/docs/user-documentation/api-reference/executors/dask). This is a parameter that you can control. For instance, you may want to define the executor to be based on Slurm using the [`SlurmExecutor`](https://docs.covalent.xyz/docs/user-documentation/api-reference/executors/slurm) to submit a job to an HPC cluster. The example below highlights how one can change the executor.

    **Setting Executors via the Lattice Object**

    If you want to use the same executor for all the `Electron` objects in a `Lattice`, you can pass the `executor` keyword argument to the `@ct.lattice` decorator, as shown below.

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc.recipes.emt.core import relax_job, static_job

    @ct.lattice(executor="local")
    def workflow(atoms):

        result1 = relax_job(atoms)
        result2 = static_job(result1)

        return result2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    **Setting Executors via the Electron Objects**

    The individual `Electron` executor options can be modified after they are imported as follows:

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc.recipes.emt.core import relax_job, static_job

    @ct.lattice
    def workflow(atoms):
        job1 = relax_job
        job1.electron_object.executor = "dask" # (1)!

        job2 = static_job
        job2.electron_object.executor = "local"

        output1 = job1(atoms)
        output2 = job2(output1)
        return output2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    1.  If you are defining your own workflow functions to use, you can also set the executor for individual `Electron` objects by passing the `executor` keyword argument to the `@ct.electron` decorator.

    **Configuring Executors**

    !!! Tip

        Refer to the [executor documentation](https://docs.covalent.xyz/docs/features/executor-plugins/exe) for instructions on how to configure Covalent for your desired machines.

    By default, the `workdir` for the `Dask` (default) and `local` executors is set to `~/.cache/covalent/workdir`. This is where any files generated at runtime will be stored. You can change both of these parameters to the directories of your choosing by editing the Covalent configuration file directly or via the `ct.set_config()` command.

    For submitting jobs to [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/) from your local machine, an example `SlurmExecutor` configuration with support for an [`sshproxy`](https://docs.nersc.gov/connect/mfa/#sshproxy)-based multi-factor authentication certificate might look like the following:

    ```python
    import covalent as ct
    n_nodes = 1
    n_cores_per_node = 48

    executor = ct.executor.SlurmExecutor(
        username="YourUserName",
        address="perlmutter-p1.nersc.gov",
        ssh_key_file="~/.ssh/nersc",
        cert_file="~/.ssh/nersc-cert.pub",
        remote_workdir="$SCRATCH/quacc",
        conda_env="quacc",
        options={
            "nodes": f"{n_nodes}",
            "qos": "debug",
            "constraint": "cpu",
            "account": "YourAccountName",
            "job-name": "quacc",
            "time": "00:10:00",
        },
        prerun_commands=[
            f"source ~/.bashrc && export QUACC_VASP_PARALLEL_CMD='srun -N {n_nodes} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores'",
        ],
        use_srun=False, # (1)!
    )
    ```

    1.  The `SlurmExecutor` must have `use_srun=False` in order for ASE-based calculators to be launched appropriately.

=== "Parsl"

    Out-of-the-box, Parsl will run on your local machine. However, in practice you will probably want to run your Parsl workflows on HPC machines.

    **Configuring Executors**

    !!! Tip

        To configure Parsl for the high-performance computing environment of your choice, refer to the executor [Configuration](https://parsl.readthedocs.io/en/stable/userguide/configuring.html) page in the Parsl documentation.

    For [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/), example [`HighThroughputExecutor`](https://parsl.readthedocs.io/en/stable/stubs/parsl.executors.HighThroughputExecutor.html#parsl.executors.HighThroughputExecutor) configurations can be found in the [NERSC Documentation](https://docs.nersc.gov/jobs/workflow/parsl/). A simple one is reproduced below that allows for job submission from the login node. This example will create a single Slurm job that will run one `PythonApp` at a time on a single node and is good for testing out some of the examples above.

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    config = Config(
        max_idletime=120, # (1)!
        executors=[
            HighThroughputExecutor(
                label="quacc_HTEX", # (2)!
                max_workers=1, # (3)!
                provider=SlurmProvider( # (4)!
                    account="MyAccountName", # (5)!
                    nodes_per_block=1, # (6)
                    scheduler_options="#SBATCH -q debug -C cpu", # (7)!
                    worker_init="source ~/.bashrc && conda activate quacc", # (8)!
                    walltime="00:10:00", # (9)!
                    cmd_timeout=120, # (10)!
                    launcher = SimpleLauncher(), # (11)!
                ),
            )
        ],
    )

    parsl.load(config)
    ```

    1. The maximum amount of time (in seconds) to allow the executor to be idle before the Slurm job is cancelled.

    2. A label for the executor instance, used during file I/O.

    3. Maximum number of workers to allow on a node.

    4. The provider to use for job submission. This can be changed to `LocalProvider()` if you wish to have the Parsl process run on a compute node rather than the login node.

    5. Your Slurm account name.

    6. The number of nodes to request per job. By default, all cores on the node will be requested (setting `cores_per_node` will override this).

    7. Any additional `#SBATCH` options can be included here.

    8. Commands to run before the job starts, typically used for activating a given Python environment.

    9. The maximum amount of time to allow the job to run in `HH:MM:SS` format.

    10. The maximum time to wait (in seconds) for the job scheduler info to be retrieved/sent.

    11. The type of Launcher to use. Note that `SimpleLauncher()` must be used instead of the commonly used `SrunLauncher()` to allow quacc subprocesses to launch their own `srun` commands.

    Unlike some other workflow engines, Parsl (by default) is built for "jobpacking" where the allocated nodes continually pull in new workers (until the walltime is reached or the parent Python process is killed). This makes it possible to request a large number of nodes that continually pull in new jobs rather than submitting a large number of small jobs to the scheduler, which can be more efficient. In other words, don't be surprised if the Slurm job continues to run even when your submitted task has completed, particularly if you are using a Jupyter Notebook or IPython kernel.

    **Scaling Up**

    Now let's consider a more realistic scenario. Suppose we want to have a single Slurm job that reserves 8 nodes, and each `PythonApp` (e.g. VASP calculation) will run on 2 nodes (let's assume each node has 48 cores total, so that's a total of 96 cores for each calculation). Parsl will act as an orchestrator in the background of one of the nodes. Our config will now look like the following.

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    n_parallel_calcs = 4 # Number of quacc calculations to run in parallel
    n_nodes_per_calc = 2 # Number of nodes to reserve for each calculation
    n_cores_per_node = 48 # Number of CPU cores per node
    vasp_parallel_cmd = f"srun -N {n_nodes_per_calc} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores"

    config = Config(
        max_idletime=300,
        executors=[
            HighThroughputExecutor(
                label="quacc_HTEX",
                max_workers=n_parallel_calcs,
                cores_per_worker=1e-6, # (1)!
                provider=SlurmProvider(
                    account="MyAccountName",
                    nodes_per_block=n_nodes_per_calc*n_parallel_calcs,
                    scheduler_options="#SBATCH -q debug -C cpu",
                    worker_init=f"source ~/.bashrc && conda activate quacc && module load vasp && export QUACC_VASP_PARALLEL_CMD={vasp_parallel_cmd}",
                    walltime="00:10:00",
                    launcher = SimpleLauncher(),
                    cmd_timeout=120,
                    init_blocks=0, # (2)!
                    min_blocks=1, # (3)!
                    max_blocks=1, # (4)!
                ),
            )
        ],
    )

    parsl.load(config)
    ```

    1. We set this to a small value so that the pilot job (e.g. the Parsl orchestrator) is allowed to be oversubscribed with scheduling processes.

    2. Sets the number of blocks (e.g. Slurm jobs) to provision during initialization of the workflow. We set this to 0 so that we only begin queuing once a workflow is submitted.

    3. Sets the minimum number of blocks (e.g. Slurm jobs) to maintain during [elastic resource management](https://parsl.readthedocs.io/en/stable/userguide/execution.html#elasticity). We set this to 1 so that the pilot job is always running. 

    4. Sets the maximum number of active blocks (e.g. Slurm jobs) during [elastic resource management](https://parsl.readthedocs.io/en/stable/userguide/execution.html#elasticity). We set this to 1 here for demonstration purposes, but it can be increased to have multiple Slurm jobpacks running simultaneously.

    !!! Tip

        Dr. Logan Ward has a nice example on YouTube describing a very similar example [here](https://youtu.be/0V4Hs4kTyJs?t=398).

=== "Prefect"

    Out-of-the-box, Prefect will run on your local machine. However, in practice you will probably want to run your Prefect workflows on HPC machines.


    **Defining Task Runners**

    !!! Tip

        Check out the [Task Runner](https://docs.prefect.io/latest/concepts/task-runners/) documentation for more information on how Prefect handles task execution.

    To modify where tasks are run, set the `task_runner` keyword argument of the corresponding `@flow` decorator. The jobs in this scenario would be submitted from a login node.

    An example is shown below for setting up a task runner compatible with the NERSC Perlmutter machine. By default, [`quacc.util.dask.make_runner`](https://quantum-accelerators.github.io/quacc/reference/quacc/util/dask.html#quacc.util.dask.make_runner) will generate a [`prefect_dask.DaskTaskRunner`](https://prefecthq.github.io/prefect-dask/task_runners/#prefect_dask.task_runners.DaskTaskRunner) composed of a [`dask_jobqueue.SLURMCluster`](https://jobqueue.dask.org/en/latest/generated/dask_jobqueue.SLURMCluster.html) object.

    ```python
    from quacc.util.dask import make_runner

    n_slurm_jobs = 1 # Number of Slurm jobs to launch in parallel
    n_nodes_per_calc = 1 # Number of nodes to reserve for each Slurm job
    n_cores_per_node = 48 # Number of CPU cores per node
    mem_per_node = "64 GB" # Total memory per node

    cluster_kwargs = {
        # Dask worker options
        "n_workers": n_slurm_jobs, # number of Slurm jobs to launch
        "cores": n_cores_per_node, # total number of cores (per Slurm job) for Dask worker
        "memory": mem_per_node, # total memory (per Slurm job) for Dask worker
        # SLURM options
        "shebang": "#!/bin/bash",
        "account": "AccountName",
        "walltime": "00:10:00", # DD:HH:SS
        "job_mem": "0", # all memory on node
        "job_script_prologue": ["source ~/.bashrc", "conda activate quacc"], # commands to run before calculation, including exports
        "job_directives_skip": ["-n", "--cpus-per-task"], # Slurm directives we can skip
        "job_extra_directives": [f"-N {n_nodes_per_calc}", "-q debug", "-C cpu"], # num. of nodes for calc (-N), queue (-q), and constraints (-c)
        "python": "python", # Python executable name
    }

    runner = make_runner(cluster_kwargs, temporary=True)
    ```

    With this instantiated cluster object, you can set the task runner of the `Flow` as follows.

    ```python
    @flow(task_runner=cluster)
    def workflow(atoms):
        ...
    ```

    Now, when the worklow is run from the login node, it will be submitted to the job scheduling system (Slurm by default), and the results will be sent back to Prefect Cloud once completed.

    !!! Tip

        Refer to the [Dask-Jobqueue Documentation](https://jobqueue.dask.org/en/latest/generated/dask_jobqueue.SLURMCluster.html) for the available `cluster_kwargs` that can be defined and how they relate to a typical job script.


    To asynchronously spawn a Slurm job that continually pulls in work for the duration of its walltime (rather than starting and terminating over the lifetime of the associated `Flow`), you can instead use the `make_runner` command without a `temporary` keyword argument:

    ```python
    runner = make_runner(cluster_kwargs)
    ```

    This is often more efficient for running large numbers of workflows because you can request a single, large Slurm job that continually pulls in work rather than submitting a large number of small jobs to the scheduler.

    Additionally, you can have the generated Dask cluster adaptively scale based on the amount of work available by setting `adapt_kwargs` as follows:

    ```python
    runner = make_runner(cluster_kwargs, adapt_kwargs={"minimum": 1, "maximum": 5})
    ```

    This will ensure that at least one Slurm job is always running, but the number of jobs will scale up to 5 if there is enough work available.

    **Troubleshooting**

    If you are having trouble figuring out the right `cluster_kwargs` to use, the best option is to have the generated job script printed to the screen. This can be done as follows.

    ```python
    from dask_jobqueue import SLURMCluster
    from quacc.util.dask import _make_cluster

    cluster = make_cluster(SLURMCluster, cluster_kwargs, verbose=True)
    ```

    Note, however, that a Slurm job will be immediately submitted, so you will probably want to `scancel` it as you debug your job script.

    **Executor Configuration File**

    Speaking of configurations, if you use mostly the same HPC settings for your calculations, it can be annoying to define a large dictionary in every workflow you run. Instead, you can define a configuration file at `~/.config/dask/jobqueue.yaml` as described in the [dask-jobqueue documentation](https://jobqueue.dask.org/en/latest/configuration-setup.html#managing-configuration-files) that can be used to define default values common to your HPC setup.

    **Using a Prefect Agent**

    So far, we have dispatched calculations immediately upon calling them. However, in practice, it is often more useful to have a [Prefect agent](https://docs.prefect.io/concepts/work-pools/#agent-overview) running in the background that will continually poll for work to submit to the task runner. This allows you to submit only a subset of workflows at a time, and the agent will automatically submit more jobs as the resources become available. You will want to run Prefect workflows with an agent on the computing environment where you wish to submit jobs, specifically on a perpetual resource like a login node or dedicated workflow node.

    To create a new work pool, run the following from the command line:

    ```bash
    prefect worker start
    ```

    It is often easiest to have this process be run in the background. Once the work pool is created, an agent can be initiated as follows:

    ```bash
    prefect agent start -p 'work-pool-name'
    ```

=== "Jobflow"

    Out-of-the box, Jobflow can be used to run on your local machine. You will, however, need a "manager" to run your workflows on HPC machines. The currently recommended manager for Jobflow is FireWorks, which is described here.

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

    The `my_qadapter.yaml` file you made in the [installation instructions](../install/install.md) specifies how FireWorks will submit jobs added to your launch pad. Additional details can be found in the [Jobflow Documentation](https://materialsproject.github.io/jobflow/tutorials/8-fireworks.html#setting-where-jobs-are-dispatched) for how to dynamically set where and how Jobflow `Job` and `Flow` objects can be dispatched.

    **Dispatching Calculations**

    With a workflow added to your launch pad, on the desired machine of choice, you can run `qlaunch rapidfire --nlaunches <N>` (where `<N>` is the number of jobs to submit) in the command line to submit your workflows to the job scheduler. Running `qlaunch rapidfire -m <N>` will ensure that `<N>` jobs are always in the queue or running. To modify the order in which jobs are run, a priority can be set via `lpad set_priority <priority> -i <FWID>` where `<priority>` is a number.

    By default, `qlaunch` will launch compute jobs that each poll for a single FireWork to run. This means that more Slurm jobs may be submitted than there are jobs to run. To modify the behavior of `qlaunch` to only submit a Slurm job for each "READY" FireWork in the launchpad, use the `-r` ("reserved") flag.

    **Monitoring the Launchpad**

    The easiest way to monitor the state of your launched FireWorks and workflows is through the GUI, which can be viewed with `lpad webgui`. To get the status of running fireworks from the command line, you can run `lpad get_fws -s RUNNING`. Other statuses can also be provided as well as individual FireWorks IDs.

    To rerun a specific FireWork, one can use the `rerun_fws` command like so: `lpad rerun_fws -i <FWID>` where `<FWID>` is the FireWork ID. Similarly, one can rerun all fizzled jobs via `lpad rerun_fws -s FIZZLED`. More complicated Mongo-style queries can also be carried out. Cancelling a workflow can be done with `lpad delete_wflows -i <FWID>`.

    Refer to the `lpad -h` help menu for more details.

    **Continuous Job Submission**

    To ensure that jobs are continually submitted to the queue you can use `tmux` to preserve the job submission process even when the SSH session is terminated. For example, running `tmux new -s launcher` will create a new `tmux` session named `launcher`. To exit the `tmux` session while still preserving any running tasks on the login node, press `ctrl+b` followed by `d`. To re-enter the tmux session, run `tmux attach -t launcher`. Additional `tmux` commands can be found on the [tmux cheatsheet](https://tmuxcheatsheet.com/).
