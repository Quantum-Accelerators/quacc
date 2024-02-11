# Worked Examples on HPC

In this section, we provide a few examples going through the entire process to deploy recipes remotely on HPC machines that use a job scheduler. The precise configuration details will depend on your given compute setup. Nonetheless, we have provided examples here for [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/) that you can build from.

!!! Tip "First-Time Deployment"

    Before deploying remote calculations for the first time, do `quacc set WORKFLOW_ENGINE null` on the remote machine and run your recipe as a standard Python script (e.g. by submitting it as a job to the scheduler). This preliminary test will help you identify potential issues early on. When you're done, you can re-set the `WORKFLOW_ENGINE` variable and continue with deployment via a workflow manager.

## Pre-Requisites

If you haven't done so already:

=== "Covalent"

    On both the local and remote machines:

    ```bash
    pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[covalent]
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

    Also make sure to connect to Prefect Cloud if you are not self-hosting:

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

    1. FireWorks has its own mechanisms for task isolation, which we will rely on instead.

## Example

When deploying calculations for the first time, it's important to start simple, which is why you should try to run a sample EMT workflow first.

=== "Covalent"

    Run the following code on the local machine:

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job

    username = "MyUserName"
    account = "MyAccountName"

    executor = ct.executor.HPCExecutor(
        username=username,
        address="perlmutter-p1.nersc.gov",
        ssh_key_file="~/.ssh/nersc",
        cert_file="~/.ssh/nersc-cert.pub",
        instance="slurm",
        resource_spec_kwargs={
            "node_count": 1,
            "processes_per_node": 1,
        },
        job_attributes_kwargs={
            "duration": 10,
            "project_name": account,
            "custom_attributes": {"slurm.constraint": "cpu", "slurm.qos": "debug"},
        },
        remote_conda_env="quacc",
        remote_workdir="$SCRATCH/quacc",
        create_unique_workdir=True,
        cleanup=False,
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

    ??? Tip "Debugging"

        The most common cause of issues is related to the job scheduler details (i.e. the `resource_spec_kwargs` and the `job_attributes_kwargs`). If your job fails on the remote machine, check the files left behind in the working directory as well as the `~/.psij` directory for a history and various log files associated with your attempted job submissions.

=== "Dask"

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

        The `dask-jobqueue` documentation has a [very helpful section](https://jobqueue.dask.org/en/latest/advanced-tips-and-tricks.html#how-to-handle-job-queueing-system-walltime-killing-workers) on how to ensure that workflows run to completion despite the finite walltime on a job scheduler system. Namely, you should use the `--lifetime` (and, potentially, the `--lifetime-stagger`) option alongside an adaptive Dask cluster to ensure that the cluster can always spin up new workers as-needed.

=== "Parsl"

    **Starting Small**

    We will start with a minimal example to get started. No concurrency will be achieved here, but it will demonstrate the basic setup.

    From an interactive resource like a Jupyter Notebook or IPython kernel on the remote machine:

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"

    config = Config(
        executors=[
            HighThroughputExecutor(
                label="quacc_parsl",
                provider=SlurmProvider(
                    account=account,
                    qos="debug",
                    constraint="cpu",
                    worker_init="source ~/.bashrc && conda activate quacc",
                    walltime="00:10:00",
                    nodes_per_block=1,
                    launcher=SimpleLauncher(),
                    cmd_timeout=120,
                ),
            )
        ],
    )

    parsl.load(config)
    ```

    ```python
    from ase.build import bulk
    from quacc.recipes.emt.core import relax_job, static_job


    def workflow(atoms):
        relax_output = relax_job(atoms)
        return static_job(relax_output["atoms"])


    atoms = bulk("Cu")
    future = workflow(atoms)
    result = future.result()
    print(result)
    ```

    **Scaling Up**

    Now it's time to scale things up and show off Parsl's true power. Let's run a TBLite relaxation and frequency calculation for 162 molecules in the so-called "g2" collection of small, neutral molecules.

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

    Now we define the workflow:

    ```python
    from quacc.recipes.tblite.core import relax_job, freq_job


    def workflow(atoms):
        relax_output = relax_job(atoms)
        return freq_job(relax_output["atoms"], energy=relax_output["results"]["energy"])
    ```

    We now loop over all molecules in the "g2" collection and apply our workflow.

    ```python
    from ase.collections import g2

    futures = []
    for name in g2.names:
        atoms = g2[name]
        future = workflow(atoms)  #  (1)!
        futures.append(future)
    ```

    1. This is where the calculations are asynchronously launched.

    We monitor the progress of our calculations and print a few summary values.

    ```python
    from tqdm import tqdm
    from concurrent.futures import as_completed

    for future in tqdm(as_completed(futures), total=len(futures)):
        task_doc = future.result()
        print(
            task_doc["formula_pretty"],
            task_doc["results"]["gibbs_energy"],
            task_doc["dir_name"],
        )
    ```

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

    Then run the following on the remote machine:

    ```bash
    qlaunch rapidfire -m 1
    ```

For an example involving a code with more complex settings, refer to the [VASP example](../advanced/vasp_hpc.md).
