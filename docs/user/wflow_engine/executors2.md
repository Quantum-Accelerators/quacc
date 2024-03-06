# Worked Examples on HPC

In this section, we provide a few examples going through the entire process to deploy recipes remotely on HPC machines that use a job scheduler. The precise configuration details will depend on your given compute setup. Nonetheless, we have provided examples here for [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/) that you can build from.

!!! Tip "First-Time Deployment"

    Before deploying remote calculations for the first time, do `quacc set WORKFLOW_ENGINE None` on the remote machine and run your recipe as a standard Python script (e.g. by submitting it as a job to the scheduler). This preliminary test will help you identify potential issues early on. When you're done, you can re-set the `WORKFLOW_ENGINE` variable and continue with deployment via a workflow manager.

## Pre-Requisites

If you haven't done so already:

=== "Covalent"

    On both the local and remote machines:

    ```bash
    pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[covalent] covalent-slurm-plugin
    quacc set WORKFLOW_ENGINE covalent && quacc set CREATE_UNIQUE_DIR false  # (1)!
    ```

    1. Many of the Covalent executors have their own mechanism for task isolation, which we will use instead.

    On the local machine:

    ```bash
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
    pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[dask]
    quacc set WORKFLOW_ENGINE dask
    ```

=== "Parsl"

    On the remote machine:

    ```bash
    pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[parsl]
    quacc set WORKFLOW_ENGINE parsl
    ```

=== "Prefect"

    On the remote machine:

    ```bash
    pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[prefect]
    quacc set WORKFLOW_ENGINE prefect
    ```

    Also make sure to connect to Prefect Cloud on the remote machine if you are not self-hosting:

    ```bash
    prefect cloud login
    ```

=== "Jobflow"

    On both the local and remote machines:

    ```bash
    pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[jobflow]
    quacc set WORKFLOW_ENGINE jobflow && quacc set CREATE_UNIQUE_DIR false  # (1)!
    ```

    1. FireWorks and Jobflow have their own mechanisms for task isolation, which we will rely on instead.

## Example

When deploying calculations for the first time, it's important to start simple, which is why you should try to run a sample EMT workflow first.

=== "Covalent"

    Run the following code on the local machine. Calculations will be dispatched to the remote machine automatically via SSH.

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job

    n_nodes = 1
    n_cores_per_node = 128

    executor = ct.executor.SlurmExecutor(
        username="YourUserName",
        address="perlmutter-p1.nersc.gov",
        ssh_key_file="/home/username/.ssh/nersc",
        cert_file="/home/username/.ssh/nersc-cert.pub",  # (1)!
        conda_env="quacc",
        options={
            "nodes": f"{n_nodes}",
            "qos": "debug",
            "constraint": "cpu",
            "account": "YourAccountName",
            "job-name": "quacc",
            "time": "00:10:00",
        },
        remote_workdir="/path/to/workdir",  # (2)!
        create_unique_workdir=True,  # (3)!
        use_srun=False,  # (4)!
        cleanup=False,  # (5)!
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

    4.  The `SlurmExecutor` must have `use_srun=False` in order for ASE-based calculators to be launched appropriately.

    5. For debugging purposes, it can be useful to keep all the temporary files. Once you're confident things work, you can omit the `cleanup` keyword argument.

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

        Then you can use the `HPCExecutor` as follows:

        ```python
        n_nodes = 1  # Number of nodes to reserve for each calculation
        n_cores_per_node = 128  # Number of CPU cores per node

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

        4.  If you use this keyword argument, there is no need to explicitly specify the `RESULTS_DIR` quacc setting.

        5. This will tell Covalent to make a unique working directory for each job. This should be used in place of the `CREATE_UNIQUE_DIR` quacc setting, which seeks to do largely the same thing.

        6. For debugging purposes, it can be useful to keep all the temporary files. Once you're confident things work, you can omit the `cleanup` keyword argument.

=== "Dask"

    From an interactive resource like a Jupyter Notebook or IPython kernel on the login node of the remote machine, run the following to instantiate a Dask [`SLURMCluster`](https://jobqueue.dask.org/en/latest/generated/dask_jobqueue.SLURMCluster.html):

    ```python
    from dask.distributed import Client
    from dask_jobqueue import SLURMCluster

    n_slurm_jobs = 1
    n_nodes_per_calc = 1
    n_cores_per_node = 128
    mem_per_node = "64 GB"
    account = "MyAccountName"

    cluster_kwargs = {
        # Dask worker options
        "n_workers": n_slurm_jobs,
        "cores": n_cores_per_node,
        "memory": mem_per_node,
        # SLURM options
        "shebang": "#!/bin/bash",
        "account": account,
        "walltime": "00:10:00",
        "job_mem": "0",
        "job_script_prologue": [
            "source ~/.bashrc",
            "conda activate quacc",
        ],
        "job_directives_skip": ["-n", "--cpus-per-task"],
        "job_extra_directives": [f"-N {n_nodes_per_calc}", "-q debug", "-C cpu"],
        "python": "python",
    }

    cluster = SLURMCluster(**cluster_kwargs)
    client = Client(cluster)
    ```

    Then run the following code:

    ```python
    from ase.build import bulk
    from quacc.recipes.emt.core import relax_job, static_job


    def workflow(atoms):
        relax_output = relax_job(atoms)
        return static_job(relax_output["atoms"])


    atoms = bulk("Cu")
    delayed = workflow(atoms)
    result = client.submit(delayed).result()
    print(result)
    ```

    !!! Tip "Handling the Walltime Killing Workers"

        The `dask-jobqueue` documentation has a [helpful section](https://jobqueue.dask.org/en/latest/advanced-tips-and-tricks.html#how-to-handle-job-queueing-system-walltime-killing-workers) on how to ensure that workflows run to completion despite the finite walltime on a job scheduler system. Namely, you should use the `--lifetime` (and, potentially, the `--lifetime-stagger`) option alongside an adaptive Dask cluster to ensure that the cluster can always spin up new workers as-needed.

=== "Parsl"

    Here, we describe several representative [`HighThroughputExecutor`](https://parsl.readthedocs.io/en/stable/stubs/parsl.executors.HighThroughputExecutor.html#parsl.executors.HighThroughputExecutor) configurations that will orchestrate jobs from the login node of NERSC's Perlmutter machine. There is no one-size-fits-all approach, so you will need to adjust the configuration to suit your specific needs.

    **Concurrent Non-MPI Jobs**

    Let's imagine a scenario where we want to concurrently run a large number of single-core compute tasks. A sample configuration for this purpose is shown below. Here, we are requesting a single Slurm allocation with 2 nodes and each job is running on one core of that allocation.

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"

    concurrent_jobs = 128
    cores_per_node = 64
    min_slurm_allocations = 0
    max_slurm_allocations = 1

    config = Config(
        strategy="htex_auto_scale",  # (1)!
        executors=[
            HighThroughputExecutor(
                label="quacc_parsl",  # (2)!
                max_workers=cores_per_node,  # (3)!
                provider=SlurmProvider(
                    account=account,
                    qos="debug",
                    constraint="cpu",
                    worker_init=f"source ~/.bashrc && conda activate quacc",  # (4)!
                    walltime="00:10:00",  # (5)!
                    nodes_per_block=concurrent_jobs // cores_per_node,  # (6)!
                    init_blocks=0,  # (7)!
                    min_blocks=min_slurm_allocations,  # (8)!
                    max_blocks=max_slurm_allocations,  # (9)!
                    launcher=SimpleLauncher(),  # (10)!
                    cmd_timeout=60,  # (11)!
                ),
            )
        ],
    )

    parsl.load(config)
    ```

    1. Unique to the `HighThroughputExecutor`, this `strategy` will automatically scale the number of active blocks (i.e. Slurm allocations) up or down based on the number of tasks remaining. We set `max_blocks=1` here so it can't scale up beyond 1 Slurm job, but it can scale down from 1 to 0 since `min_blocks=0`. By setting `init_blocks=0`, no Slurm allocation will be requested until tasks are launched.

    2. This is just an arbitrary label for file I/O.

    3. If you are running a single-core `Job`, this value will be the number of physical cores per node (this example). If you are running a `Job` that is calling MPI and uses multiple cores per node, this will be the number of concurrent tasks (see next example).

    4. Any commands to run before carrying out any of the Parsl tasks. This is useful for setting environment variables, activating a given Conda environment, and loading modules.

    5. The walltime for each block (i.e. Slurm allocation).

    6. The number of nodes that each block (i.e. Slurm allocation) should allocate.

    7. Sets the number of blocks (e.g. Slurm allocations) to provision during initialization of the workflow. We set this to a value of 0 so that there isn't a running Slurm job before any tasks have been submitted to Parsl.

    8. Sets the minimum number of blocks (e.g. Slurm allocations) to maintain during [elastic resource management](https://parsl.readthedocs.io/en/stable/userguide/execution.html#elasticity). We set this to 0 so that Slurm jobs aren't running when there are no remaining tasks.

    9. Sets the maximum number of active blocks (e.g. Slurm allocations) during [elastic resource management](https://parsl.readthedocs.io/en/stable/userguide/execution.html#elasticity). We set this to 1 here, but it can be increased to have multiple Slurm jobs running simultaneously. Raising `max_blocks` to a larger value will allow the "htex_auto_scale" strategy to upscale resources as needed.

    10. The type of Launcher to use. `SimpleLauncher()` must be used instead of the commonly used `SrunLauncher()` to allow quacc subprocesses to launch their own `srun` commands.

    11. The maximum time to wait (in seconds) for the job scheduler info to be retrieved/sent.

    **Concurrent MPI Jobs**

    Now let's consider a similar configuration but for tasks where the underlying executable is run via MPI, as is typically the case for most quantum chemistry codes that distribute work over multiple cores and/or nodes. The setup here is a bit different. In this example, we are requesting a single Slurm allocation with 8 nodes and each job is running on 2 nodes of that allocation.

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"

    concurrent_jobs = 4
    nodes_per_job = 2  # (1)!
    cores_per_node = 64  # (2)!
    min_slurm_allocations = 0
    max_slurm_allocations = 1

    config = Config(
        strategy="htex_auto_scale",
        executors=[
            HighThroughputExecutor(
                label="quacc_parsl",
                max_workers=concurrent_jobs,  # (3)!
                cores_per_worker=1e-6,  # (4)!
                provider=SlurmProvider(
                    account=account,
                    qos="debug",
                    constraint="cpu",
                    worker_init=f"source ~/.bashrc && conda activate quacc",
                    walltime="00:10:00",
                    nodes_per_block=concurrent_jobs * nodes_per_job,
                    init_blocks=0,
                    min_blocks=min_slurm_allocations,
                    max_blocks=max_slurm_allocations,
                    launcher=SimpleLauncher(),
                    cmd_timeout=60,
                ),
            )
        ],
    )

    parsl.load(config)
    ```

    1. The underlying executable should be run via `-N 2` in the `srun` command.

    2. The underlying executable should be run via `--ntasks-per-node 64` in the `srun` command. If you want the job to use half of the compute node's resources, you would set `--ntasks-per-node 32` in the MPI call but would keep `cores_per_node=64` in the Parsl configuration.

    3. Unlike the prior example, here `max_workers` is defining how many concurrent tasks to run and not how many tasks are run per node.

    4. This is recommended in the Parsl manual for jobs that run via MPI.

    **Demonstration**

    Here, we will run a TBLite relaxation and frequency calculation for 162 molecules in the so-called "g2" collection of small, neutral molecules. This example should be run from an interactive resource like a Jupyter Notebook or IPython kernel on the remote machine:

    On the remote machine, make sure to install the necessary dependencies:

    ```
    pip install quacc[tblite]
    ```

    First we initialize a Parsl configuration. For this example, we will request one Slurm job (block), which will run single-core compute tasks over two nodes.

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"

    concurrent_jobs = 128
    cores_per_node = 64
    min_slurm_allocations = 0
    max_slurm_allocations = 1

    config = Config(
        strategy="htex_auto_scale",
        executors=[
            HighThroughputExecutor(
                label="quacc_parsl",
                max_workers=cores_per_node,
                provider=SlurmProvider(
                    account=account,
                    qos="debug",
                    constraint="cpu",
                    worker_init="source ~/.bashrc && conda activate quacc",
                    walltime="00:10:00",
                    nodes_per_block=concurrent_jobs // cores_per_node,
                    init_blocks=0,
                    min_blocks=min_slurm_allocations,
                    max_blocks=max_slurm_allocations,
                    launcher=SimpleLauncher(),
                    cmd_timeout=60,
                ),
            )
        ],
    )
    parsl.load(config)
    ```

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
        task_doc = future.result()
        print(
            task_doc["formula_pretty"],
            task_doc["results"]["gibbs_energy"],
            task_doc["dir_name"],
        )
    ```

    !!! Note "Practical Deployment"

        For debugging purposes or when running only a small numbers of jobs, it is simple enough to run the Parsl process from an interactive Jupyter Notebook or IPython kernel on the remote machine. However, for practical deployment and to ensure jobs are continually submitted to the queue even when the SSH session is terminated, you can run the Parsl orchestration process on a login node and maintain its state via a program like `tmux` or `screen`.

        For example, running `tmux new -s launcher` will create a new `tmux` session named `launcher`. To exit the `tmux` session while still preserving any running tasks on the login node, press `ctrl+b` followed by `d`. To re-enter the tmux session, run `tmux attach -t launcher`. Additional `tmux` commands can be found on the [tmux cheatsheet](https://tmuxcheatsheet.com/).

=== "Prefect"

    From an interactive resource like a Jupyter Notebook or IPython kernel on the login node of the remote machine, run the following to instantiate a Dask [`SLURMCluster`](https://jobqueue.dask.org/en/latest/generated/dask_jobqueue.SLURMCluster.html):

    ```python
    from dask.distributed import Client
    from dask_jobqueue import SLURMCluster

    n_slurm_jobs = 1
    n_nodes_per_calc = 1
    n_cores_per_node = 48
    mem_per_node = "64 GB"

    cluster_kwargs = {
        # Dask worker options
        "n_workers": n_slurm_jobs,
        "cores": n_cores_per_node,
        "memory": mem_per_node,
        # SLURM options
        "shebang": "#!/bin/bash",
        "account": "MyAccountName",  # (1)!
        "walltime": "00:10:00",
        "job_mem": "0",
        "job_script_prologue": [
            "source ~/.bashrc",
            "conda activate quacc",
        ],
        "job_directives_skip": ["-n", "--cpus-per-task"],
        "job_extra_directives": [f"-N {n_nodes_per_calc}", "-q debug", "-C cpu"],
        "python": "python",
    }

    cluster = SLURMCluster(**cluster_kwargs)
    client = Client(cluster)
    ```

    1. Make sure to replace this with the account name to charge.

    Then run the following code:

    ```python
    from ase.build import bulk
    from prefect_dask import DaskTaskRunner
    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job


    @flow(task_runner=DaskTaskRunner(address=client.scheduler.address))
    def workflow(atoms):
        relax_output = relax_job(atoms)
        return static_job(relax_output["atoms"])


    atoms = bulk("Cu")
    future = workflow(atoms)
    print(future.result())
    ```

    If you are connecting to an existing Dask cluster and want to ensure it is not killed when the walltime is reached, refer to the [corresponding section](https://jobqueue.dask.org/en/latest/advanced-tips-and-tricks.html#how-to-handle-job-queueing-system-walltime-killing-workers) in the `dask-jobqueue` manual. If preferred, it is also possible to instantiate a [one-time, temporary](https://prefecthq.github.io/prefect-dask/usage_guide/#using-a-temporary-cluster) Dask cluster via the `DaskTaskRunner` rather than connecting to an existing Dask cluster. This is the more conventional job scheduling approach, where each workflow will run on its own job allocation.

=== "Jobflow"

    From the login node of the remote machine, run the following:

    ```python
    import jobflow as jf
    from ase.build import bulk
    from fireworks import LaunchPad
    from jobflow.managers.fireworks import flow_to_workflow
    from quacc.recipes.emt.core import relax_job, static_job

    atoms = bulk("Cu")
    job1 = relax_job(atoms)
    job2 = static_job(job1.output["atoms"])
    flow = jf.Flow([job1, job2])

    wf = flow_to_workflow(flow)
    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf)
    ```

    **Dispatching Calculations**

    With a workflow added to your launch pad, on the desired machine of choice, you can run `qlaunch rapidfire --nlaunches <N>` (where `<N>` is the number of jobs to submit) in the command line to submit your workflows to the job scheduler. Running `qlaunch rapidfire -m <N>` will ensure that `<N>` jobs are always in the queue or running. To modify the order in which jobs are run, a priority can be set via `lpad set_priority <priority> -i <FWID>` where `<priority>` is a number.

    By default, `qlaunch` will launch compute jobs that each poll for a single FireWork to run. This means that more Slurm jobs may be submitted than there are jobs to run. To modify the behavior of `qlaunch` to only submit a Slurm job for each "READY" FireWork in the launchpad, use the `-r` ("reserved") flag.

    **Monitoring the Launchpad**

    The easiest way to monitor the state of your launched FireWorks and workflows is through the GUI, which can be viewed with `lpad webgui`. To get the status of running fireworks from the command line, you can run `lpad get_fws -s RUNNING`. Other statuses can also be provided as well as individual FireWorks IDs.

    To rerun a specific FireWork, one can use the `rerun_fws` command like so: `lpad rerun_fws -i <FWID>` where `<FWID>` is the FireWork ID. Similarly, one can rerun all fizzled jobs via `lpad rerun_fws -s FIZZLED`. More complicated Mongo-style queries can also be carried out. Cancelling a workflow can be done with `lpad delete_wflows -i <FWID>`.

    Refer to the `lpad -h` help menu for more details.

    **Continuous Job Submission**

    To ensure that jobs are continually submitted to the queue, you can use `tmux` to preserve the job submission process even when the SSH session is terminated. For example, running `tmux new -s launcher` will create a new `tmux` session named `launcher`. To exit the `tmux` session while still preserving any running tasks on the login node, press `ctrl+b` followed by `d`. To re-enter the tmux session, run `tmux attach -t launcher`. Additional `tmux` commands can be found on the [tmux cheatsheet](https://tmuxcheatsheet.com/).

For an example involving a code with more complex settings, refer to the [VASP example](../advanced/vasp_hpc.md).
