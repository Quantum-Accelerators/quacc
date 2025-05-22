# Deploying Calculations

In the previous examples, we have been running calculations on the local machine. However, in practice, you will probably want to run your calculations on one or more HPC machines. This section will describe how to set up your workflows to run on HPC machines using your desired workflow engine to scale up your calculations.

!!! Note "A Note on Terminology"

    Throughout this section, when we refer to a "job", we are specifically referring to a function decorated by `#!Python @job` (i.e. not necessarily a Slurm job). A node refers to the compute node on the HPC machine, which consists of many individual compute cores. An allocation refers to one or more nodes that have been reserved via the job scheduler for a given amount of time.

## Background Information

=== "Covalent"

    Refer to the [executor plugin documentation](https://docs.covalent.xyz/docs/plugin) for instructions on how to install and use the relevant plugins that allow Covalent to submit jobs on your desired machines. Most users of quacc will probably want to use the [`SlurmExecutor`](https://github.com/AgnostiqHQ/covalent-slurm-plugin), which is a plugin for Covalent that supports Slurm job scheduling system. If you are using a job scheduler environment but find the `covalent-slurm-plugin` does not suit your needs, you may wish to consider the [`covalent-hpc-plugin`](https://github.com/Quantum-Accelerators/covalent-hpc-plugin).

=== "Dask"

    A Dask cluster can be set up to be used with a queueing system like that found on most HPC machines. This is most easily done via [Dask Jobqueue](https://jobqueue.dask.org/en/latest/index.html). Example configurations for various queuing systems can be found in the ["Example Deployments"](https://jobqueue.dask.org/en/latest/examples.html) section of the Dask Jobqueue documentation.

=== "Parsl"

    Out-of-the-box, Parsl will run on your local machine. However, in practice you will probably want to run your Parsl workflows on HPC machines.

    To configure Parsl for the high-performance computing environment of your choice, refer to the [executor configuration page in the Parsl documentation](https://parsl.readthedocs.io/en/stable/userguide/configuring.html) for many examples. Additional details can be found in the ["Execution" section](https://parsl.readthedocs.io/en/stable/userguide/execution.html) of the Parsl documentation. Most users of quacc will probably want to the [`HighThroughputExecutor`](https://parsl.readthedocs.io/en/stable/stubs/parsl.executors.HighThroughputExecutor.html#parsl.executors.HighThroughputExecutor).

    !!! Note "Pilot Jobs"

        Unlike most other workflow engines, Parsl is built for the [pilot job model](https://en.wikipedia.org/wiki/Pilot_job) where the allocated nodes continually pull in new jobs to run. This makes it possible to avoid submitting a large number of small jobs to the scheduler, which can be inefficient from a queuing perspective.

    !!! Tip "Globus Compute"

        If you want to run Parsl workflows on remote, distributed computing resources, check out [Globus Compute](https://github.com/globus/globus-compute) and the corresponding [tutorial](https://globus-compute.readthedocs.io/en/latest/tutorial.html#running-parsl-workflows).

=== "Prefect"

    To scale up calculations, read about the concept of a Prefect [task runner](https://docs.prefect.io/latest/concepts/task-runners/). By default, `quacc` automatically submits all `#!Python @job`-decorated functions to the specified task runner and so concurrency is achieved by default.

    To use Prefect in a job scheduler environment, you can create a [`DaskTaskRunner`](https://prefecthq.github.io/prefect-dask/usage_guide/) that can be used in conjunction with [dask-jobqueue](https://jobqueue.dask.org/en/latest). Example configurations for various queuing systems can be found in the ["Example Deployments"](https://jobqueue.dask.org/en/latest/examples.html) section of the `dask-jobqueue` documentation.

=== "Jobflow"

    Out-of-the-box, Jobflow can be used to run on your local machine. You will, however, need a "manager" to run your workflows on HPC machines. The recommended option is [jobflow-remote](https://github.com/Matgenix/jobflow-remote), but you can also use [Fireworks](https://github.com/materialsproject/fireworks).

    === "Jobflow Remote"

        **Setting Up Your Jobflow Remote Cnnfiguration**

        When you previously [set up Jobflow or FireWorks](../../install/wflow_engines.md), you created a YAML file with configuration details. It's now time to revisit that file and adjust the `pre_run` command with any modules or environment variables necessary for your calculations to run. Additionally, you will probably want to update the `nodes`, `walltime`, and related settings for your scheduler.

    === "Fireworks"

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

## Worked Examples

In this section, we go through the entire process to deploy recipes remotely on HPC machines that use a job scheduler. The precise configuration details will depend on your given compute setup. Nonetheless, we have provided examples here for [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/) that you can build from.

!!! Tip "First-Time Deployment"

    Before deploying remote calculations for the first time, do `quacc set WORKFLOW_ENGINE None` on the remote machine and run your recipe as a standard Python script (e.g. by submitting it as a job to the scheduler). This preliminary test will help you identify potential issues that are independent of the workflow manager.

### Pre-Requisites

If you haven't done so already:

=== "Covalent"

    On both the local and remote machines:

    ```bash
    pip install 'quacc[covalent]'
    quacc set WORKFLOW_ENGINE covalent && quacc set CREATE_UNIQUE_DIR false  # (1)!
    ```

    1. Many of the Covalent executors have their own mechanism for task isolation, which we will use instead.

    On the local machine:

    ```bash
    pip install covalent-slurm-plugin
    covalent start
    ```

    ??? Note "For NERSC Users"

        If using Perlmutter at NERSC, modify your `~/.bashrc` on the remote machine as follows since only the `$SCRATCH` directory supports file locking mechanisms:

        ```bash title="~/.bashrc"
        export COVALENT_CONFIG_DIR="$SCRATCH/.config/covalent"
        ```

=== "Dask"

    On the remote machine:

    ```bash
    pip install 'quacc[dask]'
    quacc set WORKFLOW_ENGINE dask
    ```

    ??? Note "For NERSC Users"

        If using Perlmutter at NERSC, Dask should be run from the `$SCRATCH` directory. This is because the `$SCRATCH` directory is the only directory that supports file locking mechanisms, which Dask relies on.

=== "Parsl"

    On the remote machine:

    ```bash
    pip install 'quacc[parsl]'
    quacc set WORKFLOW_ENGINE parsl
    ```

=== "Prefect"

    On the remote machine:

    ```bash
    pip install 'quacc[prefect]'
    quacc set WORKFLOW_ENGINE prefect
    ```

    Also make sure to connect to Prefect Cloud on the remote machine if you are not self-hosting:

    ```bash
    prefect cloud login
    ```

    ??? Note "For NERSC Users"

        If using Perlmutter at NERSC with the Dask backend for Prefect, your calculations should be run from the `$SCRATCH` directory. This is because the `$SCRATCH` directory is the only directory that supports file locking mechanisms, which Dask relies on.

=== "Jobflow"

    On both the local and remote machines:

    === "Jobflow Remote"

        ```bash
        pip install 'quacc[jobflow]'
        ```

    === "Fireworks"

        ```bash
        pip install 'quacc[jobflow]' fireworks
        ```

### Concurrent Non-MPI Jobs

=== "Covalent"

    Run the following code on the local machine. Calculations will be dispatched to the remote machine automatically via SSH.

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job

    nodes = 1

    executor = ct.executor.SlurmExecutor(
        username="YourUserName",
        address="perlmutter-p1.nersc.gov",
        ssh_key_file="~/.ssh/nersc",
        cert_file="~/.ssh/nersc-cert.pub",  # (1)!
        conda_env="quacc",
        options={
            "nodes": nodes,
            "qos": "debug",
            "constraint": "cpu",
            "account": "YourAccountName",
            "job-name": "quacc",
            "time": "00:10:00",
        },
        remote_workdir="/path/to/workdir",  # (2)!
        create_unique_workdir=True,  # (3)!
        cleanup=False,  # (4)!
    )


    @flow(executor=executor)
    def workflow(atoms):
        relax_output = relax_job(atoms)
        return static_job(relax_output["atoms"])


    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    1.  This a certificate file used to validate your SSH credentials. This is often not needed but is required at NERSC facilities due to the use of [`sshproxy`](https://docs.nersc.gov/connect/mfa/#sshproxy)-based multi-factor authentication.

    2.  If you use this keyword argument, there is no need to explicitly specify the `RESULTS_DIR` quacc setting.

    3.  This will tell Covalent to make a unique working directory for each job. This should be used in place of the `CREATE_UNIQUE_DIR` quacc setting, which seeks to do largely the same thing.

    4. For debugging purposes, it can be useful to keep all the temporary files. Once you're confident things work, you can omit the `cleanup` keyword argument.

    ??? Note "An Alternate Approach: The `HPCExecutor`"

        If you are using the `HPCExecutor`, the process is similar.

        On the local machine, run:

        ```bash
        pip install covalent-hpc-plugin
        ```

        On the remote machine, run:

        ```bash
        pip install psij-python
        ```

        Then you can use the `HPCExecutor` as follows from the local machine:

        ```python
        import covalent as ct

        nodes = 1
        cores_per_node = 128

        executor = ct.executor.HPCExecutor(
            # SSH credentials
            username="YourUserName",
            address="perlmutter-p1.nersc.gov",
            ssh_key_file="~/.ssh/nersc",
            cert_file="~/.ssh/nersc-cert.pub",  # (1)!
            # PSI/J parameters
            instance="slurm",
            resource_spec_kwargs={
                "node_count": nodes,
                "processes_per_node": cores_per_node,
            },  # (2)!
            job_attributes_kwargs={
                "duration": 10,  # minutes
                "account": "YourAccountName",
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

        4.  If you use this keyword argument, there is no need to explicitly specify the `RESULTS_DIR` quacc setting.

        5. This will tell Covalent to make a unique working directory for each job. This should be used in place of the `CREATE_UNIQUE_DIR` quacc setting, which seeks to do largely the same thing.

        6. For debugging purposes, it can be useful to keep all the temporary files. Once you're confident things work, you can omit the `cleanup` keyword argument.

=== "Dask"

    Here, we will run single-core TBLite relaxation and frequency calculations for 162 molecules in the so-called "g2" collection of small, neutral molecules.

    On the remote machine, first make sure to install the necessary dependencies:

    ```
    pip install 'quacc'[tblite]'
    ```

    From an interactive resource like a Jupyter Notebook or IPython kernel on the login node of the remote machine, run the following to instantiate a Dask [`SLURMCluster`](https://jobqueue.dask.org/en/latest/generated/dask_jobqueue.SLURMCluster.html) that will request two Slurm allocations of one node each, with each `#!Python @job` running on one core of the allocation:

    ```python
    from dask.distributed import Client
    from dask_jobqueue import SLURMCluster

    account = "MyAccountName"

    cores_per_node = 128
    slurm_jobs = 2

    env_vars = "export OMP_NUM_THREADS=1,1"  # (1)!

    cluster = SLURMCluster(
        cores=cores_per_node,  # (2)!
        memory="64 GB",
        shebang="#!/bin/bash",
        account=account,
        walltime="00:10:00",
        job_mem="0",
        job_script_prologue=[
            "source ~/.bashrc",
            env_vars,
        ],
        job_directives_skip=["-n", "--cpus-per-task"],  # (3)!
        job_extra_directives=["-q debug", "-C cpu"],  # (4)!
    )
    print(cluster.job_script())
    cluster.scale(jobs=slurm_jobs)
    client = Client(cluster)
    ```

    1. Since we are running single-core jobs, we need to set the `OMP_NUM_THREADS` environment variable to "1,1" according to the [TBLite documentation](https://tblite.readthedocs.io/en/latest/tutorial/parallel.html#running-tblite-in-parallel).

    2. It is also worthwhile to play around with the `processes` argument, which defaults to sqrt(cores).

    3. We have skipped the `-n` and `--cpus-per-task` Slurm directives because we will be using the full node here.

    4. We have set the queue to submit to and the constraint to use.

    Then run the following code:

    ```python
    from ase.collections import g2
    from quacc.recipes.tblite.core import relax_job, freq_job


    def workflow(atoms):
        relax_output = relax_job(atoms)
        return freq_job(relax_output["atoms"], energy=relax_output["results"]["energy"])


    futures = []
    for name in g2.names:
        atoms = g2[name]
        future = client.compute(workflow(atoms))
        futures.append(future)

    results = client.gather(futures)
    for result in results:
        print(
            result["formula_pretty"],
            result["results"]["gibbs_energy"],
            result["dir_name"],
        )
    ```

=== "Parsl"

    Here, we describe representative [`HighThroughputExecutor`](https://parsl.readthedocs.io/en/stable/stubs/parsl.executors.HighThroughputExecutor.html#parsl.executors.HighThroughputExecutor) configurations that will orchestrate jobs from the login node of NERSC's Perlmutter machine. There is no one-size-fits-all approach, so you will need to adjust the configuration to suit your specific needs.

    Here, we will run single-core TBLite relaxation and frequency calculations for 162 molecules in the so-called "g2" collection of small, neutral molecules. This example should be run from an interactive resource like a Jupyter Notebook or IPython kernel on the remote machine.

    On the remote machine, make sure to install the necessary dependencies:

    ```
    pip install 'quacc[tblite]'
    ```

    Now we will request a single Slurm allocation with 2 nodes, and each compute job will run on one core of that allocation.

    ```python
    import parsl
    from parsl.config import Config
    from parsl.dataflow.dependency_resolvers import DEEP_DEPENDENCY_RESOLVER
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SrunLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"

    cores_per_job = 1
    cores_per_node = 128
    nodes_per_allocation = 2
    min_allocations = 0
    max_allocations = 1

    env_vars = f"export OMP_NUM_THREADS={cores_per_job},1"  # (1)!

    config = Config(
        dependency_resolver=DEEP_DEPENDENCY_RESOLVER,  # (2)!
        strategy="htex_auto_scale",  # (3)!
        executors=[
            HighThroughputExecutor(
                label="quacc_parsl",  # (4)!
                max_workers_per_node=cores_per_node,  # (5)!
                cores_per_worker=cores_per_job,  # (6)!
                provider=SlurmProvider(
                    account=account,
                    qos="debug",
                    constraint="cpu",
                    worker_init=f"source ~/.bashrc && conda activate cms && {env_vars}",  # (7)!
                    walltime="00:10:00",  # (8)!
                    nodes_per_block=nodes_per_allocation,  # (9)!
                    init_blocks=0,  # (10)!
                    min_blocks=min_allocations,  # (11)!
                    max_blocks=max_allocations,  # (12)!
                    launcher=SrunLauncher(),  # (13)!
                    cmd_timeout=60,  # (14)!
                ),
            )
        ],
        initialize_logging=False,  # (15)!
    )

    parsl.load(config)
    ```

    1. Since we are running single-core jobs, we need to set the `OMP_NUM_THREADS` environment variable to "1,1" according to the [TBLite documentation](https://tblite.readthedocs.io/en/latest/tutorial/parallel.html#running-tblite-in-parallel).

    2. This is an opt-in feature needed when using Parsl to ensure all features in quacc are supported.

    3. Unique to the `HighThroughputExecutor`, this `strategy` will automatically scale the number of active blocks (i.e. Slurm allocations) up or down based on the number of jobs remaining. We set `max_blocks=1` here so it can't scale up beyond 1 Slurm job, but it can scale down from 1 to 0 since `min_blocks=0`. By setting `init_blocks=0`, no Slurm allocation will be requested until jobs are launched.

    4. This is just an arbitrary label for file I/O.

    5. The maximum number of running jobs per node. If you are running a non-MPI job, this value will generally be the number of physical cores per node (this example). Perlmutter has 128 physical CPU cores, so we have set a value of 128 here.

    6. The number of cores per job. We are running single-core jobs in this example.

    7. Any commands to run before carrying out any of the Parsl jobs. This is useful for setting environment variables, activating a given Conda environment, and loading modules.

    8. The walltime for each block (i.e. Slurm allocation).

    9. The number of nodes that each block (i.e. Slurm allocation) should allocate.

    10. Sets the number of blocks (e.g. Slurm allocations) to provision during initialization of the workflow. We set this to a value of 0 so that there isn't a running Slurm job before any jobs have been submitted to Parsl.

    11. Sets the minimum number of blocks (e.g. Slurm allocations) to maintain during [elastic resource management](https://parsl.readthedocs.io/en/stable/userguide/execution.html#elasticity). We set this to 0 so that Slurm jobs aren't running when there are no remaining jobs.

    12. Sets the maximum number of active blocks (e.g. Slurm allocations) during [elastic resource management](https://parsl.readthedocs.io/en/stable/userguide/execution.html#elasticity). We set this to 1 here, but it can be increased to have multiple Slurm jobs running simultaneously. Raising `max_blocks` to a larger value will allow the "htex_auto_scale" strategy to upscale resources as needed.

    13. The type of Launcher to use. `SrunLauncher()` will distribute jobs across the cores and nodes of the Slurm allocation. It should not be used for `PythonApp`s that themselves call MPI, which should use `SimpleLauncher()` instead.

    14. The maximum time to wait (in seconds) for the job scheduler info to be retrieved/sent.

    15. This will tidy up the Parsl logging to match the same log level as in quacc (`INFO` by default).

    Now we define the workflow, apply it to all molecules in the "g2" collection, and monitor the progress of our calculations.

    ```python
    from ase.collections import g2
    from concurrent.futures import as_completed
    from tqdm import tqdm
    from quacc.recipes.tblite.core import relax_job, freq_job


    def workflow(atoms):
        relax_output = relax_job(atoms)
        return freq_job(relax_output["atoms"], energy=relax_output["results"]["energy"])


    futures = []
    for name in g2.names:
        atoms = g2[name]
        future = workflow(atoms)  #  (1)!
        futures.append(future)


    for future in tqdm(as_completed(futures), total=len(futures)):
        result = future.result()
        print(
            result["formula_pretty"],
            result["results"]["gibbs_energy"],
            result["dir_name"],
        )
    ```

    1. This is when the `PythonApp`s will be dispatched.

    ??? Tip "Tips for Practical Deployment"

        For debugging purposes or when running only a small numbers of jobs, it is simple enough to run the Parsl process from an interactive Jupyter Notebook or IPython kernel on the remote machine. However, for practical deployment and to ensure jobs are continually submitted to the queue even when the SSH session is terminated, you can run the Parsl orchestration process on a login node and maintain its state via a program like `tmux` or `screen`.

        For example, running `tmux new -s launcher` will create a new `tmux` session named `launcher`. To exit the `tmux` session while still preserving any running jobs on the login node, press `ctrl+b` followed by `d`. To re-enter the tmux session, run `tmux attach -t launcher`. Additional `tmux` commands can be found on the [tmux cheatsheet](https://tmuxcheatsheet.com/).

=== "Prefect"

    Here, we will run single-core TBLite relaxation and frequency calculations for 20 molecules from the so-called "g2" collection of small, neutral molecules. Note that the full "g2" collection has 162 molecules, but running all of them in quick succession will risk surpassing the free-tier rate limit of Prefect Cloud.

    On the remote machine, first make sure to install the necessary dependencies:

    ```
    pip install 'quacc[tblite]'
    ```

    From an interactive resource like a Jupyter Notebook or IPython kernel on the login node of the remote machine, you will need to instantiate a Dask [`SLURMCluster`](https://jobqueue.dask.org/en/latest/generated/dask_jobqueue.SLURMCluster.html) that will request one Slurm allocation for one node, with a `#!Python @job` running on each core of the allocation:

    ```python
    from dask.distributed import Client
    from dask_jobqueue import SLURMCluster

    account = "MyAccountName"

    slurm_jobs = 1
    cores_per_node = 128

    env_vars = "export OMP_NUM_THREADS=1,1"  # (1)!

    cluster_kwargs = {
        "cores": cores_per_node,  # (2)!
        "memory": "64 GB",
        "shebang": "#!/bin/bash",
        "account": account,
        "walltime": "00:10:00",
        "job_mem": "0",
        "job_script_prologue": [
            "source ~/.bashrc",
            env_vars,
        ],
        "job_directives_skip": ["-n", "--cpus-per-task"],  # (3)!
        "job_extra_directives": ["-q debug", "-C cpu"],  # (4)!
    }
    cluster = SLURMCluster(**cluster_kwargs)
    print(cluster.job_script())
    cluster.scale(jobs=slurm_jobs)
    client = Client(cluster)
    ```

    1. Since we are running single-core jobs, we need to set the `OMP_NUM_THREADS` environment variable to "1,1" according to the [TBLite documentation](https://tblite.readthedocs.io/en/latest/tutorial/parallel.html#running-tblite-in-parallel).

    2. It is also worthwhile to play around with the `processes` argument, which defaults to sqrt(cores).

    3. We have skipped the `-n` and `--cpus-per-task` Slurm directives because we will be using the full node here.

    4. We have set the queue to submit to and the constraint to use.

    Now we define our workflow to dispatch, attaching it to the premade Dask cluster:

    ```python
    from prefect_dask import DaskTaskRunner
    from quacc import flow
    from quacc.recipes.tblite.core import freq_job, relax_job


    @flow(task_runner=DaskTaskRunner(address=client.scheduler.address))
    def workflow(atoms_objects):
        futures = []
        for atoms in atoms_objects:
            relax_output = relax_job(atoms)
            freq_output = freq_job(
                relax_output["atoms"], energy=relax_output["results"]["energy"]
            )
            futures.append(freq_output)

        return futures
    ```

    Finally, we dispatch the workflow and fetch the results:

    ```python
    from ase.collections import g2

    atoms_objects = [g2[name] for name in g2.names[:20]]
    results = workflow(atoms_objects)
    ```

    !!! Tip "One-Time Dask Clusters"

        If preferred, it is also possible to instantiate a [one-time, temporary](https://prefecthq.github.io/prefect-dask/usage_guide/#using-a-temporary-cluster) Dask cluster via the `DaskTaskRunner` directly rather than connecting to an existing Dask cluster, as described in the [Task Runner documentation](https://prefecthq.github.io/prefect-dask/task_runners/). This is the more conventional job scheduling approach, where each workflow will run on its own Slurm allocation.

=== "Jobflow"

    === "Jobflow Remote"

        From the login node of the remote machine, run the following:

        ```python
        from ase.build import bulk
        from jobflow import Flow
        from jobflow_remote import submit_flow
        from quacc.recipes.emt.core import relax_job, static_job

        atoms_list = [bulk("Si"), bulk("Al")]
        for atoms in atoms_list:
            job1 = relax_job(atoms, relax_cell=True)
            job2 = static_job(job1.output["atoms"])
            flow = Flow([job1, job2])

            submit_flow(flow, worker="basic_python")
        ```

    === "Fireworks"

        From the login node of the remote machine, run the following:

        ```python
        import jobflow as jf
        from ase.build import bulk
        from fireworks import LaunchPad
        from jobflow.managers.fireworks import flow_to_workflow
        from quacc.recipes.emt.core import relax_job, static_job

        lpad = LaunchPad.auto_load()

        atoms_list = [bulk("Si"), bulk("Al")]
        for atoms in atoms_list:
            job1 = relax_job(atoms, relax_cell=True)
            job2 = static_job(job1.output["atoms"])
            flow = Flow([job1, job2])

            wf = flow_to_workflow(flow)
            lpad.add_wf(wf)
        ```

        **Dispatching Calculations**

        With a workflow added to your launch pad, on the desired machine of choice, you can run `qlaunch rapidfire --nlaunches <N>` (where `<N>` is the number of jobs to submit) in the command line to submit your workflows to the job scheduler. Running `qlaunch rapidfire -m <N>` will ensure that `<N>` jobs are always in the queue or running. To modify the order in which jobs are run, a priority can be set via `lpad set_priority <priority> -i <FWID>` where `<priority>` is a number.

        By default, `qlaunch` will launch compute jobs that each poll for a single FireWork to run. This means that more Slurm jobs may be submitted than there are jobs to run. To modify the behavior of `qlaunch` to only submit a Slurm job for each "READY" FireWork in the launchpad, use the `-r` ("reserved") flag.

        **Monitoring the Launchpad**

        The easiest way to monitor the state of your launched FireWorks and workflows is through the GUI, which can be viewed with `lpad webgui`. To get the status of running fireworks from the command line, you can run `lpad get_fws -s RUNNING`. Other statuses can also be provided as well as individual FireWorks IDs.

        To rerun a specific FireWork, one can use the `rerun_fws` command like so: `lpad rerun_fws -i <FWID>` where `<FWID>` is the FireWork ID. Similarly, one can rerun all fizzled jobs via `lpad rerun_fws -s FIZZLED`. More complicated Mongo-style queries can also be carried out. Cancelling a workflow can be done with `lpad delete_wflows -i <FWID>`. Refer to the `lpad -h` help menu for more details.

        **Setting Where Jobs are Dispatched**

        The `my_qadapter.yaml` file you made in the [installation instructions](../../install/install.md) specifies how FireWorks will submit jobs added to your launch pad. Additional details can be found in the [Jobflow Documentation](https://materialsproject.github.io/jobflow/tutorials/8-fireworks.html#setting-where-jobs-are-dispatched) for how to dynamically set where and how Jobflow `Job` and `Flow` objects can be dispatched.

        ??? Tip "Continuous Job Submission"

            To ensure that jobs are continually submitted to the queue, you can use `tmux` to preserve the job submission process even when the SSH session is terminated. For example, running `tmux new -s launcher` will create a new `tmux` session named `launcher`. To exit the `tmux` session while still preserving any running jobs on the login node, press `ctrl+b` followed by `d`. To re-enter the tmux session, run `tmux attach -t launcher`. Additional `tmux` commands can be found on the [tmux cheatsheet](https://tmuxcheatsheet.com/).

### Concurrent MPI Jobs

Here we will run a sample VASP recipe that will highlight the use of a more complicated MPI-based configuration. This example can only be run if you are a licensed VASP user, but the same fundamental principles apply to many other DFT codes with recipes in quacc.

First, prepare your `QUACC_VASP_PP_PATH` environment variable in the `~/.bashrc` of your remote machine as described in the [Calculator Setup guide](../../install/codes.md). When you're done, follow the steps below.

=== "Covalent"

    Run the following code on the local machine:

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.vasp.core import relax_job, static_job

    nodes = 1
    cores_per_node = 128
    vasp_parallel_cmd = (
        f"srun -N {nodes} --ntasks-per-node={cores_per_node} --cpu_bind=cores"
    )

    executor = ct.executor.SlurmExecutor(
        username="YourUserName",
        address="perlmutter-p1.nersc.gov",
        ssh_key_file="~/.ssh/nersc",
        cert_file="~/.ssh/nersc-cert.pub",
        conda_env="quacc",
        options={
            "nodes": nodes,
            "qos": "debug",
            "constraint": "cpu",
            "account": "YourAccountName",
            "job-name": "quacc",
            "time": "00:30:00",
        },
        remote_workdir="/path/to/workdir",
        create_unique_workdir=True,
        use_srun=False,  # (1)!
        prerun_commands=[
            "module load vasp/6.4.1-cpu",
            f"export QUACC_VASP_PARALLEL_CMD='{vasp_parallel_cmd}'",
        ],
    )


    @flow(executor=executor)
    def workflow(atoms):
        relax_output = relax_job(atoms, kpts=[3, 3, 3])
        return static_job(relax_output["atoms"], kpts=[3, 3, 3])


    atoms = bulk("C")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    1. This is set to `False` because we are using the `srun` command in the `QUACC_VASP_PARALLEL_CMD` environment variable.

    ??? Note "An Alternate Approach: The `HPCExecutor`"

        A similar configuration can be used for the `HPCExecutor`:

        ```python
        import covalent as ct

        nodes = 1
        cores_per_node = 128
        vasp_parallel_cmd = (
            f"srun -N {nodes} --ntasks-per-node={cores_per_node} --cpu_bind=cores"
        )

        executor = ct.executor.HPCExecutor(
            username="YourUserName",
            address="perlmutter-p1.nersc.gov",
            ssh_key_file="~/.ssh/nersc",
            cert_file="~/.ssh/nersc-cert.pub",
            instance="slurm",
            resource_spec_kwargs={
                "node_count": nodes,
                "processes_per_node": cores_per_node,
            },
            job_attributes_kwargs={
                "duration": 30,
                "account": "YourAccountName",
                "custom_attributes": {"slurm.constraint": "cpu", "slurm.qos": "debug"},
            },
            pre_launch_cmds=["module load vasp/6.4.1-cpu"],
            environment={
                "QUACC_VASP_PARALLEL_CMD": vasp_parallel_cmd,
            },
            remote_conda_env="quacc",
            remote_workdir="$SCRATCH/quacc",
            create_unique_workdir=True,
            cleanup=False,
        )
        ```

=== "Dask"

    ```python
    from dask.distributed import Client
    from dask_jobqueue import SLURMCluster

    account = "MyAccountName"

    slurm_jobs = 2
    nodes_per_calc = 1
    cores_per_node = 128

    vasp_parallel_cmd = (
        f"srun -N {nodes_per_calc} --ntasks-per-node={cores_per_node} --cpu_bind=cores"
    )

    cluster = SLURMCluster(
        cores=1,  # (1)!
        memory="64 GB",
        shebang="#!/bin/bash",
        account=account,
        walltime="00:10:00",
        job_mem="0",
        job_script_prologue=[
            "source ~/.bashrc",
            "module load vasp/6.4.1-cpu",
            f"export QUACC_VASP_PARALLEL_CMD='{vasp_parallel_cmd}'",
        ],
        job_directives_skip=["-n", "--cpus-per-task"],
        job_extra_directives=[f"-N {nodes_per_calc}", "-q debug", "-C cpu"],
    )
    print(cluster.job_script())
    cluster.scale(jobs=slurm_jobs)
    client = Client(cluster)
    ```

    1. Since we only want to run one VASP job per allocation, we set this to 1. VASP will still use multiple cores via the `srun` command. The same could be achieved by setting `processes=1`.

    Now we launch the VASP jobs:

    ```python
    from ase.build import bulk
    from quacc.recipes.vasp.core import relax_job, static_job


    def workflow(atoms):
        relax_output = relax_job(atoms, kpts=[3, 3, 3])
        return static_job(relax_output["atoms"], kpts=[3, 3, 3])


    future1 = client.compute(workflow(bulk("C")))
    future2 = client.compute(workflow(bulk("Cu")))
    result = client.gather([future1, future2])
    print(result)
    ```

    !!! Warning "Limitations"

        Unfortunately, `dask-jobqueue` is somewhat limited in terms of its flexibility. Most notably, there is no mechanism to distribute `#!Python @job`s over multiple nodes on a single Slurm allocation. Users interested in such functionality should consider using other tools in the Dask suite or other workflow engines (e.g. Parsl).

=== "Parsl"

    Now let's consider a similar configuration but for jobs where the underlying executable is run via MPI, as is typically the case for most quantum chemistry codes that distribute work over multiple cores and/or nodes. The setup here is a bit different. In this example, we are requesting a single Slurm allocation with 2 nodes (containing 128 physical CPU cores per node), and each compute job is running on 1 node of that allocation.

    From an interactive resource like a Jupyter Notebook or IPython kernel from the login node on the remote machine:

    ```python
    import parsl
    from parsl.config import Config
    from parsl.dataflow.dependency_resolvers import DEEP_DEPENDENCY_RESOLVER
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"

    nodes_per_job = 1
    cores_per_node = 128
    nodes_per_allocation = 2
    vasp_parallel_cmd = (
        f"srun -N {nodes_per_job} --ntasks-per-node={cores_per_node} --cpu_bind=cores"
    )

    config = Config(
        dependency_resolver=DEEP_DEPENDENCY_RESOLVER,
        strategy="htex_auto_scale",
        executors=[
            HighThroughputExecutor(
                label="quacc_mpi_parsl",
                max_workers_per_node=nodes_per_allocation // nodes_per_job,  # (1)!
                cores_per_worker=1e-6,  # (2)!
                provider=SlurmProvider(
                    account=account,
                    qos="debug",
                    constraint="cpu",
                    worker_init=f"source ~/.bashrc && conda activate cms && module load vasp/6.4.1-cpu && export QUACC_VASP_PARALLEL_CMD='{vasp_parallel_cmd}'",
                    walltime="00:10:00",
                    nodes_per_block=nodes_per_allocation,
                    launcher=SimpleLauncher(),  # (3)!
                    cmd_timeout=60,
                ),
            )
        ],
        initialize_logging=False,
    )

    parsl.load(config)
    ```

    1. Unlike the prior example, here `max_workers_per_node` is defining the maximum number of concurrent MPI jobs to run per allocation.

    2. This is recommended in the Parsl manual for jobs that spawn MPI processes.

    3. The `SimpleLauncher` should be used in place of the `SrunLauncher` for `PythonApp`s that themselves call MPI.

    Now we launch the VASP jobs:

    ```python
    from ase.build import bulk
    from quacc.recipes.vasp.core import relax_job, static_job


    def workflow(atoms):
        relax_output = relax_job(atoms, kpts=[3, 3, 3])
        return static_job(relax_output["atoms"], kpts=[3, 3, 3])


    future1 = workflow(bulk("C"))
    future2 = workflow(bulk("Cu"))
    print(future1.result(), future2.result())
    ```

=== "Prefect"

    ```python
    from dask.distributed import Client
    from dask_jobqueue import SLURMCluster

    account = "MyAccountName"

    slurm_jobs = 2
    nodes_per_calc = 1
    cores_per_node = 128

    vasp_parallel_cmd = (
        f"srun -N {nodes_per_calc} --ntasks-per-node={cores_per_node} --cpu_bind=cores"
    )

    cluster_kwargs = {
        "cores": 1,  # (1)!
        "memory": "64 GB",
        "shebang": "#!/bin/bash",
        "account": account,
        "walltime": "00:10:00",
        "job_mem": "0",
        "job_script_prologue": [
            "source ~/.bashrc",
            "module load vasp/6.4.1-cpu",
            f"export QUACC_VASP_PARALLEL_CMD='{vasp_parallel_cmd}'",
        ],
        "job_directives_skip": ["-n", "--cpus-per-task"],
        "job_extra_directives": [f"-N {nodes_per_calc}", "-q debug", "-C cpu"],
    }
    cluster = SLURMCluster(**cluster_kwargs)
    print(cluster.job_script())
    cluster.scale(jobs=slurm_jobs)
    client = Client(cluster)
    ```

    1. Since we only want to run one VASP job per allocation, we set this to 1. VASP will still use multiple cores via the `srun` command. The same could be achieved by setting `processes=1`.

    Now we define our workflow to dispatch, attaching it to the premade Dask cluster:

    ```python
    from prefect_dask import DaskTaskRunner
    from quacc import flow


    @flow(task_runner=DaskTaskRunner(address=client.scheduler.address))
    def workflow(list_of_atoms):
        from quacc.recipes.vasp.core import relax_job, static_job

        futures = []
        for atoms in list_of_atoms:
            relax_output = relax_job(atoms, kpts=[3, 3, 3])
            static_output = static_job(relax_output["atoms"], kpts=[3, 3, 3])
            futures.append(static_output)

        return futures
    ```

    Finally, we dispatch the workflow and fetch the results:

    ```python
    from ase.build import bulk

    list_of_atoms = [bulk("Cu"), bulk("C")]
    results = workflow(list_of_atoms)
    ```

=== "Jobflow"

    === "Jobflow Remote"

        From the login node of the remote machine, run the following:

        ```python
        from ase.build import bulk
        from jobflow import Flow
        from jobflow_remote import submit_flow
        from quacc.recipes.vasp.core import relax_job, static_job

        atoms_list = [bulk("Si"), bulk("Al")]
        for atoms in atoms_list:
            atoms.set_initial_magnetic_moments([0.0] * len(atoms))
            job1 = relax_job(atoms, relax_cell=True)
            job2 = static_job(job1.output["atoms"])
            flow = Flow([job1, job2])

            submit_flow(flow, worker="basic_vasp")
        ```

        Then monitor the progress with `jf job list`.

    === "Fireworks"

        You will need to update your `my_qadapter.yaml` file that you made when setting up FireWorks. Specifically, ensure that the following parameters are set:

        ```yaml title="my_qadapter.yaml"
        _fw_name: CommonAdapter
        _fw_q_type: SLURM
        rocket_launch: rlaunch -w </path/to/fw_config/my_fworker.yaml> singleshot
        nodes: 1
        walltime: 00:30:00
        account: MySlurmAccountName
        job_name: quacc_firework
        qos: debug
        pre_rocket: |
                    conda activate cms
                    module load vasp/6.4.1-cpu
                    export QUACC_VASP_PARALLEL_CMD="srun -N 1 --ntasks-per-node=128 --cpu_bind=cores"
                    export QUACC_WORKFLOW_ENGINE=jobflow
                    export QUACC_CREATE_UNIQUE_DIR=False
        ```

        From the login node of the remote machine, run the following:

        ```python
        import jobflow as jf
        from ase.build import bulk
        from fireworks import LaunchPad
        from jobflow.managers.fireworks import flow_to_workflow
        from quacc.recipes.vasp.core import relax_job, static_job

        lpad = LaunchPad.auto_load()

        atoms_list = [bulk("Si"), bulk("Al")]
        for atoms in atoms_list:
            atoms.set_initial_magnetic_moments([0.0] * len(atoms))
            job1 = relax_job(atoms, relax_cell=True)
            job2 = static_job(job1.output["atoms"])
            flow = Flow([job1, job2])

            wf = flow_to_workflow(flow)
            lpad.add_wf(wf)
        ```

        Then run the following on the remote machine:

        ```bash
        qlaunch rapidfire -m 1
        ```
