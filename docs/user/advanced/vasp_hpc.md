# VASP Example on HPC

As a complement to the [toy EMT example](../wflow_engine/executors2.md) described earlier in the documentation, here we will run a sample VASP recipe that will highlight the use of a more complicated configuration. It can only be run if you are a licensed VASP user, but the same fundamental principles apply to many other DFT codes with recipes in quacc.

First, prepare your `QUACC_VASP_PP_PATH` environment variable in the `~/.bashrc` of your remote machine as described in the [Calculator Setup guide](../../install/codes.md). When you're done, follow the steps below.

=== "Covalent"

    Run the following code on the local machine:

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.vasp.core import relax_job, static_job

    username = "MyUserName"
    account = "MyAccountName"
    n_nodes = 1
    n_cores_per_node = 128

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
        pre_launch_cmds=["module load vasp/6.4.1-cpu"],
        environment={
            "QUACC_VASP_PARALLEL_CMD": f"srun -N {n_nodes} --ntasks-per-node={n_cores_per_node} --cpu_bind=cores"
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


    atoms = bulk("C")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    print(result)
    ```

=== "Parsl"

    Let's consider a scenario where we want to run concurrent VASP jobs that each run on a full CPU node of 128 cores. We will run two concurrent VASP jobs in a single Slurm allocation.

    From an interactive resource like a Jupyter Notebook or IPython kernel from the login node on the remote machine:

    ```python
    import parsl
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    from parsl.launchers import SimpleLauncher
    from parsl.providers import SlurmProvider

    account = "MyAccountName"

    concurrent_jobs = 2
    nodes_per_job = 1
    cores_per_node = 128
    vasp_parallel_cmd = f"srun -N {nodes_per_job} --ntasks-per-node={cores_per_node} --cpu_bind=cores"
    min_slurm_allocations = 0
    max_slurm_allocations = 1

    config = Config(
        strategy="htex_auto_scale",
        executors=[
            HighThroughputExecutor(
                label="quacc_parsl",
                max_workers=concurrent_jobs,
                cores_per_worker=1e-6,
                provider=SlurmProvider(
                    account=account,
                    qos="debug",
                    constraint="cpu",
                    worker_init=f"source ~/.bashrc && conda activate quacc && module load vasp/6.4.1-cpu && export QUACC_VASP_PARALLEL_CMD={vasp_parallel_cmd}",
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
    export QUACC_VASP_PARALLEL_CMD="srun -N 1 --ntasks-per-node=128 --cpu_bind=cores"
    ```

    From the login node of the remote machine, then run the following:

    ```python
    import jobflow as jf
    from ase.build import bulk
    from fireworks import LaunchPad
    from jobflow.managers.fireworks import flow_to_workflow
    from quacc.recipes.vasp.core import relax_job, static_job

    atoms = bulk("C")
    job1 = relax_job(atoms, kpts=[3, 3, 3])
    job2 = static_job(job1.output["atoms"], kpts=[3, 3, 3])
    flow = jf.Flow([job1, job2])

    wf = flow_to_workflow(flow)
    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf)
    ```

    Then run the following on the remote machine:

    ```bash
    qlaunch rapidfire -m 1
    ```
