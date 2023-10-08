# Worked Examples on HPC

In this section, we provide a few examples going through the entire process to deploy recipes remotely on HPC machines that use a job scheduler. The precise configuration details will depend on your given compute setup. Nonetheless, we have provided examples here for [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/) that you can build from.

!!! Tip

    Before deploying remote calculations for the first time, do `quacc set WORKFLOW_ENGINE local` on the remote machine and run your recipe as a standard Python script (e.g. by submitting it as a job to the scheduler). This preliminary test will help you identify potential issues early on. When you're done, you can re-set the `WORKFLOW_ENGINE` variable and continue with deployment via a workflow manager.

## Pre-Requisites

Start with a clean Conda environment if you don't have one already:

```bash
conda create --name quacc python=3.10
conda activate quacc
```

Then install the necessary dependencies:

=== "Covalent ⭐"

    On both the local and remote machines:

    ```bash
    pip install --no-cache-dir https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[covalent]
    quacc set WORKFLOW_ENGINE covalent
    ```

    On the local machine:

    ```bash
    covalent start
    ```

    !!! Note

        If using Perlmutter at NERSC, modify your `~/.bashrc` on the remote machine as follows since only the `$SCRATCH` directory supports file locking mechanisms:

        ```bash title="~/.bashrc"
        export COVALENT_CONFIG_DIR="$SCRATCH/.config/covalent"
        ```

=== "Parsl ⭐"

    On both the remote machine:

    ```bash
    pip install --no-cache-dir https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[parsl]
    quacc set WORKFLOW_ENGINE parsl
    quacc set CREATE_UNIQUE_WORKDIR True
    ```

=== "Jobflow"

    On both the local and remote machines:

    ```bash
    pip install --no-cache-dir https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[jobflow]
    quacc set WORKFLOW_ENGINE jobflow
    ```

## Example 1: EMT

When deploying calculations for the first time, it's important to start simple, which is why you should try to run a sample EMT workflow first.

=== "Covalent ⭐"

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


    @flow(executor=executor, workflow_executor=executor)  # (1)!
    def workflow(atoms):
        relax_output = relax_job(atoms)
        return static_job(relax_output)


    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    1. Until [Issue 1024](https://github.com/Quantum-Accelerators/quacc/issues/1024) is resolved, you need to directly set the `workflow_executor` keyword argument in the `#!Python @flow` decorator to the same value as that used for `executor` otherwise a post-processing error will occur.

    !!! Hint

        The most common cause of issues is related to the job scheduler details (i.e. the `resource_spec_kwargs` and the `job_attributes_kwargs`). If your job fails on the remote machine, check the files left behind in the working directory as well as the `~/.psij` directory for a history and various log files associated with your attempted job submissions.

=== "Parsl ⭐"

    **Starting Small**

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
                max_workers=1,
                provider=SlurmProvider(
                    account=account,
                    scheduler_options="#SBATCH -q debug -C cpu",
                    worker_init="source ~/.bashrc && conda activate quacc",
                    walltime="00:10:00",
                    nodes_per_block=1,
                    init_blocks=0,
                    min_blocks=1,
                    max_blocks=1,
                    launcher=SimpleLauncher(),
                ),
            )
        ],
    )

    parsl.load(config)
    ```

    ```python
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.emt.core import relax_job, static_job


    @flow
    def workflow(atoms):
        relax_output = relax_job(atoms)
        return static_job(relax_output)


    atoms = bulk("Cu")
    future = workflow(atoms)
    result = future.result()
    print(result)
    ```

    **Scaling Up**

    Now it's time to scale things up and show off Parsl's true power. Let's run a TBLite relaxation and frequency calculation for 162 molecules in the so-called "g2" collection of small, neutral molecules.

    On the remote machine, make sure to run `pip install quacc[tblite]`. Then run the following example.

    First we initialize a Parsl configuration. For this example, we will request 2 Slurm jobs (blocks), each of which will run tasks over 2 nodes that will be dynamically scaled

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"

    config = Config(
        max_idletime=60,
        strategy="htex_auto_scale",
        executors=[
            HighThroughputExecutor(
                label="quacc_parsl",
                provider=SlurmProvider(
                    account=account,
                    scheduler_options="#SBATCH -q debug -C cpu",
                    worker_init="source ~/.bashrc && conda activate quacc",
                    walltime="00:10:00",
                    nodes_per_block=2,
                    init_blocks=0,
                    min_blocks=0,
                    max_blocks=2,
                    launcher=SimpleLauncher(),
                    cmd_timeout=120,
                ),
            )
        ],
    )
    parsl.load(config)
    ```

    Now we define the workflow:

    ```python
    from ase.build import bulk
    from quacc.recipes.tblite.core import relax_job, freq_job
    from quacc import job

    def workflow(atoms):
        relax_output = relax_job(atoms)
        return freq_job(relax_output)
    ```

    We now loop over all molecules in the "g2" collection and apply our workflow.

    ```python
    from ase.build import molecule
    from ase.collections import g2

    futures = []
    for name in g2.names:
        atoms = molecule(name)
        future = workflow(atoms)
        futures.append(future)
    ```

    We monitor the progress of our calculations and print a few summary values.

    ```python
    import time
    from tqdm import tqdm
    from concurrent.futures import as_completed

    for future in tqdm(as_completed(futures), total=len(futures)):
        task_doc = future.result()
        print(task_doc["formula_pretty"], task_doc["results"]["gibbs_energy"], task_doc["dir_name"])
    ```

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
    job2 = static_job(job1.output)
    flow = jf.Flow([job1, job2])

    wf = flow_to_workflow(flow)
    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf)
    ```

    Then run the following on the remote machine:

    ```bash
    qlaunch rapidfire -m 1
    ```

## Example 2: VASP

In this example, we will run a sample VASP recipe that will highlight the use of a more complicated configuration.

First, prepare your `VASP_PP_PATH` environment variable in the `~/.bashrc` of your remote machine as described in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#pseudopotentials). When you're done, follow the steps below.

=== "Covalent ⭐"

    Run the following code on the local machine:

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.vasp.core import relax_job, static_job

    username = "MyUserName"
    account = "MyAccountName"
    n_nodes = 1
    n_cores_per_node = 48

    executor = ct.executor.HPCExecutor(
        username=username,
        address="perlmutter-p1.nersc.gov",
        ssh_key_file="~/.ssh/nersc",
        cert_file="~/.ssh/nersc-cert.pub",
        instance="slurm",
        resource_spec_kwargs={
            "node_count": n_nodes,
            "processes_per_node": n_cores_per_node,
        },
        job_attributes_kwargs={
            "duration": 30,
            "project_name": account,
            "custom_attributes": {"slurm.constraint": "cpu", "slurm.qos": "debug"},
        },
        pre_launch_cmds=[
            f"export QUACC_VASP_PARALLEL_CMD='srun -N {n_nodes} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores'",
            "module load vasp/6.4.1-cpu",
        ],  # (1)!
        remote_conda_env="quacc",
        remote_workdir="$SCRATCH/quacc",
        create_unique_workdir=True,
        cleanup=False,
    )


    @flow(executor=executor, workflow_executor=executor)  # (2)!
    def workflow(atoms):
        relax_output = relax_job(atoms)
        return static_job(relax_output)


    atoms = bulk("C")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    1. Until [this issue](https://github.com/ExaWorks/psij-python/issues/423) is resolved, environment variables should be specified in `pre_launch_cmds`. Once it's resolved, it can be specified as `#!Python environment={"QUACC_VASP_PARALLEL_CMD": f"srun -N {n_nodes} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores"}` instead.

    2. Until [Issue 1024](https://github.com/Quantum-Accelerators/quacc/issues/1024) is resolved, you need to directly set the `workflow_executor` keyword argument in the `#!Python @flow` decorator to the same value as that used for `executor` otherwise a post-processing error will occur.

=== "Parsl ⭐"

    From an interactive resource like a Jupyter Notebook or IPython kernel on the remote machine:

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"
    max_slurm_jobs = 1
    n_calcs_per_job = 1
    n_nodes_per_calc = 1
    n_cores_per_node = 48

    config = Config(
        strategy="htex_auto_scale",
        executors=[
            HighThroughputExecutor(
                label="quacc_parsl",
                max_workers=n_parallel_calcs,
                cores_per_worker=1e-6,
                provider=SlurmProvider(
                    account=account,
                    scheduler_options="#SBATCH -q debug -C cpu",
                    worker_init=f"source ~/.bashrc && conda activate quacc && module load vasp/6.4.1-cpu && export QUACC_VASP_PARALLEL_CMD='srun -N {n_nodes_per_calc} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores'",
                    walltime="00:10:00",
                    nodes_per_block=n_nodes_per_calc * n_calcs_per_job,
                    init_blocks=0,
                    min_blocks=0,
                    max_blocks=max_slurm_jobs,
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
    from quacc import flow
    from quacc.recipes.vasp.core import relax_job, static_job


    @flow
    def workflow(atoms):
        relax_output = relax_job(atoms, calc_swaps={"kpts": [3, 3, 3]})
        return static_job(relax_output, calc_swaps={"kpts": [3, 3, 3]})


    atoms = bulk("C")
    future = workflow(atoms)
    result = future.result()
    print(result)
    ```

=== "Jobflow"

    You will need to update your `my_qadapter.yaml` file that you made when setting up FireWorks. Specifically, ensure that the following parameters are set:

    ```yaml title="my_qadapter.yaml"
    _fw_name: CommonAdapter
    _fw_q_type: SLURM
    rocket_launch: rlaunch -w /path/to/fw_config/my_fworker.yaml singleshot
    nodes: 1
    walltime: 00:30:00
    account: MyAccountName
    job_name: quacc_firework
    qos: debug
    pre_rocket: |
    module load vasp/6.4.1-cpu
    export QUACC_VASP_PARALLEL_CMD="srun -N 1 --ntasks-per-node=48 --cpu_bind=cores"
    ```

    From the login node of the remote machine, then run the following:

    ```python
    import jobflow as jf
    from ase.build import bulk
    from fireworks import LaunchPad
    from jobflow.managers.fireworks import flow_to_workflow
    from quacc.recipes.vasp.core import relax_job, static_job

    atoms = bulk("C")
    job1 = relax_job(atoms, calc_swaps={"kpts": [3, 3, 3]})
    job2 = static_job(job1.output, calc_swaps={"kpts": [3, 3, 3]})
    flow = jf.Flow([job1, job2])

    wf = flow_to_workflow(flow)
    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf)
    ```

    Then run the following on the remote machine:

    ```bash
    qlaunch rapidfire -m 1
    ```
