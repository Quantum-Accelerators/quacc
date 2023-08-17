# Deploying Calculations

In the previous examples, we have been running calculations on our local machine. However, in practice, you will probably want to run your calculations on one or more HPC machines. This section will describe how to set up your workflows to run on HPC machines using your desired workflow engine to scale up your calculations.

=== "Covalent"

    By default, Covalent will run all `Electron` tasks on your local machine using the [`DaskExecutor`](https://docs.covalent.xyz/docs/user-documentation/api-reference/executors/dask). This is a parameter that you can control. For instance, Covalent offers many [plugin executors](https://docs.covalent.xyz/docs/features/executor-plugins/exe) that can be used to interface with a wide range of HPC, cloud, and quantum devices.

    **Setting Executors via the Lattice Object**

    If you want to use the same executor for all the `Electron` objects in a `Lattice`, you can pass the `executor` keyword argument to the `@ct.lattice` decorator, as shown below.

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc.recipes.emt.core import relax_job, static_job

    @ct.lattice(executor="local") # (1)!
    def workflow(atoms):

        result1 = relax_job(atoms)
        result2 = static_job(result1)

        return result2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    1. This was merely for demonstration purposes. There is never really a need to use the "local" executor since the "dask" executor runs locally and is faster.

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
        job2.electron_object.executor = "local" # (2)!

        output1 = job1(atoms)
        output2 = job2(output1)
        return output2

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    1. If you are defining your own workflow functions to use, you can also set the executor for individual `Electron` objects by passing the `executor` keyword argument to the `@ct.electron` decorator.

    2. This was merely for demonstration purposes. There is never really a need to use the "local" executor since the "dask" executor runs locally and is faster.

    **Configuring Executors**

    !!! Tip

        Refer to the [executor documentation](https://docs.covalent.xyz/docs/features/executor-plugins/exe) for instructions on how to install and use the relevant plugins that allow Covalent to submit jobs on your desired machines.

    Most users of quacc will probably want to use the [`HPCExecutor`](https://github.com/Quantum-Accelerators/covalent-hpc-plugin), which is a plugin for Covalent that supports Slurm, PBS, LSF, Flux, and more. For submitting jobs to [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/) from your local machine, an example `HPCExecutor` configuration with support for an [`sshproxy`](https://docs.nersc.gov/connect/mfa/#sshproxy)-based multi-factor authentication certificate might look like the following:

    ```python
    executor = ct.executor.HPCExecutor(
        # SSH credentials
        username="YourUserName",
        address="perlmutter-p1.nersc.gov",
        ssh_key_file="~/.ssh/nersc",
        cert_file="~/.ssh/nersc-cert.pub", # (1)!
        # PSI/J parameters
        instance="slurm",
        resource_spec_kwargs={
            "nodes": n_nodes,
            "processes_per_node": n_cores_per_node,
        }, # (2)!
        job_attributes_kwargs={
            "duration": 10, # minutes
            "project_name": "YourAccountName",
            "custom_attributes": {"slurm.constraint": "cpu", "slurm.qos": "debug"},
        }, # (3)!
        environment={"QUACC_VASP_PARALLEL_CMD": vasp_parallel_cmd},
        # Pre-/post-launch commands
        pre_launch_cmds=["module load vasp"],
        # Remote Python env parameters
        remote_conda_env="quacc",
        # Covalent parameters
        remote_workdir="$SCRATCH/quacc",
        create_unique_workdir=True, # (4)!
    )
    ```

    1. This a certificate file used to validate your SSH credentials. This is often not needed but is required at NERSC facilities.

    2. These are the resource specifications for the compute job, which are keyword arguments passed to PSI/J's [`ResourceSpecV1` class](https://exaworks.org/psij-python/docs/v/0.9.0/.generated/psij.html#psij.resource_spec.ResourceSpecV1).

    3. These are the job attributes that the job scheduler needs, which are keyword arguments passed to PSI/J's [`JobAttributes` class](https://exaworks.org/psij-python/docs/v/0.9.0/.generated/psij.html#psij.JobAttributes).

    4. You generally want each quacc job to be run in its own unique working directory to ensure files don't overwrite one another, so  `create_unique_workdir` should be set to `True`.

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
                    account="MyAccountName",
                    nodes_per_block=1, # (5)!
                    scheduler_options="#SBATCH -q debug -C cpu", # (6)!
                    worker_init="source ~/.bashrc && conda activate quacc",
                    walltime="00:10:00",
                    cmd_timeout=120, # (7)!
                    launcher = SimpleLauncher(), # (8)!
                ),
            )
        ],
    )

    parsl.load(config)
    ```

    1. The maximum amount of time (in seconds) to allow the executor to be idle before the Slurm job is cancelled.

    2. A label for the executor instance, used during file I/O.

    3. Maximum number of workers to allow on a node.

    4. The provider to use for job submission. This can be changed to `LocalProvider()` if you wish to have the Parsl process run on a login node rather than a compute node.

    5. The number of nodes to request per job. By default, all cores on the node will be requested (setting `cores_per_node` will override this).

    6. Any additional `#SBATCH` options not captured elsewhere can be included here.

    7. The maximum time to wait (in seconds) for the job scheduler info to be retrieved/sent.

    8. The type of Launcher to use. Note that `SimpleLauncher()` must be used instead of the commonly used `SrunLauncher()` to allow quacc subprocesses to launch their own `srun` commands.

    Unlike some other workflow engines, Parsl (by default) is built for "jobpacking" (also known as the pilot job model) where the allocated nodes continually pull in new workers (until the walltime is reached or the parent Python process is killed). This makes it possible to request a large number of nodes that continually pull in new jobs rather than submitting a large number of small jobs to the scheduler, which can be more efficient. In other words, don't be surprised if the Slurm job continues to run even when your submitted task has completed, particularly if you are using a Jupyter Notebook or IPython kernel.

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
    vasp_parallel_cmd = (
        f"srun -N {n_nodes} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores'"
    )

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

    An example is shown below for setting up a task runner compatible with the NERSC Perlmutter machine. By default, [`quacc.util.wflows.make_runner`](https://quantum-accelerators.github.io/quacc/reference/quacc/util/dask.html#quacc.util.wflows.make_runner) will generate a [`prefect_dask.DaskTaskRunner`](https://prefecthq.github.io/prefect-dask/task_runners/#prefect_dask.task_runners.DaskTaskRunner) composed of a [`dask_jobqueue.SLURMCluster`](https://jobqueue.dask.org/en/latest/generated/dask_jobqueue.SLURMCluster.html) object.

    ```python
    from quacc.util.wflows import make_runner

    n_slurm_jobs = 1 # Number of Slurm jobs to launch in parallel.
    n_nodes_per_calc = 1 # Number of nodes to reserve for each Slurm job.
    n_cores_per_node = 48 # Number of CPU cores per node.
    mem_per_node = "64 GB" # Total memory per node.
    vasp_parallel_cmd = (
        f"srun -N {n_nodes} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores'"
    )

    cluster_kwargs = {
        # Dask worker options
        "n_workers": n_slurm_jobs, # (1)!
        "cores": n_cores_per_node, # (2)!
        "memory": mem_per_node, # (3)!
        # SLURM options
        "shebang": "#!/bin/bash",
        "account": "AccountName",
        "walltime": "00:10:00",
        "job_mem": "0", # (4)!
        "job_script_prologue": [
            "source ~/.bashrc",
            "conda activate quacc",
            f"export QUACC_VASP_PARALLEL_CMD={vasp_parallel_cmd}",
        ], # (5)!
        "job_directives_skip": ["-n", "--cpus-per-task"], # (6)!
        "job_extra_directives": [f"-N {n_nodes_per_calc}", "-q debug", "-C cpu"], # (7)!
        "python": "python", # (8)!
    }

    runner = make_runner(cluster_kwargs, temporary=True)
    ```

    1. Number of Slurm jobs to launch.

    2. Total number of cores (per Slurm job) for Dask worker.

    3. Total memory (per Slurm job) for Dask worker.

    4. Request all memory on the node.

    5. Commands to run before calculation. This is a good place to include environment variable definitions and modules to load.

    6. Slurm directives that are automatically added but that we chose to skip.

    7. The number of nodes for each calculation (-N), queue name (-q), and constraint (-c). Oftentimes, the constraint flag is not needed.

    8. The Python executable name. This often does not need to be changed.

    With this instantiated cluster object, you can set the task runner of the `Flow` as follows.

    ```python
    @flow(task_runner=runner)
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
    from quacc.util.wflows import _make_cluster

    cluster = make_cluster(SLURMCluster, cluster_kwargs, verbose=True)
    ```

    Note, however, that a Slurm job will be immediately submitted, so you will probably want to `scancel` it as you debug your job script.

    **Executor Configuration File**

    Speaking of configurations, if you use mostly the same HPC settings for your calculations, it can be annoying to define a large dictionary in every workflow you run. Instead, you can define a configuration file at `~/.config/dask/jobqueue.yaml` as described in the [dask-jobqueue documentation](https://jobqueue.dask.org/en/latest/configuration-setup.html#managing-configuration-files) that can be used to define default values common to your HPC setup.

    **Using a Prefect Work Pool and Agent**

    So far, we have dispatched calculations immediately upon calling them. However, in practice, it is often more useful to have a Prefect agent running in the background that will continually poll for work to submit to the task runner. This allows you to submit only a subset of workflows at a time, and the agent will automatically submit more jobs as the resources become available. You will want to run Prefect workflows with an agent on the computing environment where you wish to submit jobs, specifically on a perpetual resource like a login node or dedicated workflow node. Refer to the ["Work Pools, Workers, and Agents"](https://docs.prefect.io/latest/concepts/work-pools/) section of the Prefect documentation for more details.
