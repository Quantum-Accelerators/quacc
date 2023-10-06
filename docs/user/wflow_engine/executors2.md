# Worked Examples

In this section, we provide a few examples going through the entire process to deploy recipes remotely for some commonly used workflow engines. The precise configuration details will depend on your given compute setup. Nonetheless, we have provided examples here for [Perlmutter at NERSC](https://docs.nersc.gov/systems/perlmutter/) that you can build from.

!!! Hint

    Before deploying remote calculations for the first time, ensure that the following can be done successfully:

    1. Run the sample recipe on your local machine, when possible.

    2. Run the same Python script on your desired computing resource (e.g. by submitting it as a job to the scheduler). Make sure that the `WORKFLOW_ENGINE` setting is set to "local" on the remote machine for this exercise.

    These preliminary tests will help you identify potential issues early on.

## Pre-Requisites

On the local and remote machines, make a clean Conda environment if you haven't already:

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

        If using Perlmutter at NERSC, modify your `~/.bashrc` on the remote machine as follows:

        ```bash title="~/.bashrc"
        export COVALENT_CONFIG_DIR="$SCRATCH/.config/covalent"
        ```

=== "Parsl ⭐"

    On both the local and remote machines:

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

    @flow(executor=executor, workflow_executor=executor) # (1)!
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

        The most common cause of issues is related to the job scheduler details (i.e. the `resource_spec_kwargs` and the `job_attributes_kwargs`). If your job fails on the remote machine, check the `~/.psij` directory for a history and various log files associated with your attempted job submissions.

=== "Parsl ⭐"

    From an interactive resource like a Jupyter Notebook or IPython kernel on the remote machine:

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"

    config = Config(
        max_idletime=300,
        executors=[
            HighThroughputExecutor(
                label="quacc_HTEX",
                max_workers=1,
                cores_per_worker=1e-6,
                provider=SlurmProvider(
                    account=account,
                    nodes_per_block=1,
                    scheduler_options="#SBATCH -q debug -C cpu",
                    worker_init=f"source ~/.bashrc && conda activate quacc",
                    walltime="00:10:00",
                    launcher=SimpleLauncher(),
                    cmd_timeout=120,
                    init_blocks=0,
                    min_blocks=1,
                    max_blocks=1,
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
    vasp_parallel_cmd = (
        f"srun -N {n_nodes} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores"
    )

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
            "duration": 10,
            "project_name": account,
            "custom_attributes": {"slurm.constraint": "cpu", "slurm.qos": "debug"},
        },
        environment={"QUACC_VASP_PARALLEL_CMD": vasp_parallel_cmd},
        pre_launch_cmds=["module load vasp"],
        remote_conda_env="quacc",
        remote_workdir="$SCRATCH/quacc",
        create_unique_workdir=True,
        cleanup=False,
    )


    @flow(executor=executor, workflow_executor=executor) # (1)!
    def workflow(atoms):
        relax_output = relax_job(atoms)
        return static_job(relax_output)

    atoms = bulk("Fe")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

    1. Until [Issue 1024](https://github.com/Quantum-Accelerators/quacc/issues/1024) is resolved, you need to directly set the `workflow_executor` keyword argument in the `#!Python @flow` decorator to the same value as that used for `executor` otherwise a post-processing error will occur.

=== "Parsl ⭐"

    From an interactive resource like a Jupyter Notebook or IPython kernel on the remote machine:

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"
    n_parallel_calcs = 1
    n_nodes_per_calc = 1
    n_cores_per_node = 48
    vasp_parallel_cmd = (
        f"srun -N {n_nodes_per_calc} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores"
    )

    config = Config(
        max_idletime=300,
        executors=[
            HighThroughputExecutor(
                label="quacc_HTEX",
                max_workers=n_parallel_calcs,
                cores_per_worker=1e-6,
                provider=SlurmProvider(
                    account=account,
                    nodes_per_block=n_nodes_per_calc * n_parallel_calcs,
                    scheduler_options="#SBATCH -q debug -C cpu",
                    worker_init=f"source ~/.bashrc && conda activate quacc && module load vasp && export QUACC_VASP_PARALLEL_CMD={vasp_parallel_cmd}",
                    walltime="00:10:00",
                    launcher=SimpleLauncher(),
                    cmd_timeout=120,
                    init_blocks=0,
                    min_blocks=1,
                    max_blocks=1,
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
        relax_output = relax_job(atoms)
        return static_job(relax_output)

    atoms = bulk("Fe")
    future = workflow(atoms)
    result = future.result()
    print(result)
    ```

=== "Jobflow"

    From the login node of the remote machine, run the following:

    ```python
    import jobflow as jf
    from ase.build import bulk
    from fireworks import LaunchPad
    from jobflow.managers.fireworks import flow_to_workflow
    from quacc.recipes.vasp.core import relax_job, static_job

    atoms = bulk("Fe")
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
