# Worked Examples

In this section, we provide a few examples going through the entire process to deploy recipes remotely for some commonly used workflow engines. The precise configuration details will depend on your given compute setup, but we have provided two examples here (one for Perlmutter at NERSC and the other for Tiger at Princeton) for you to build from.

!!! Hint

    Before deploying remote calculations for the first time, ensure that the following can be done successfully:

    1. Run the sample recipe on your local machine, where possible.

    2. Run that same Python script on your desired computing resource (e.g. by submitting it as a job to the scheduler). Make sure that the `WORKFLOW_ENGINE` setting is set to "local" on the remote machine for this exercise.

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

    === "Perlmutter"

        ```python
        import covalent as ct
        from ase.build import bulk
        from quacc import flow
        from quacc.recipes.emt.core import relax_job, static_job

        username = "MyUserName"
        address = "perlmutter-p1.nersc.gov"
        account = "MyAccountName"

        executor = ct.executor.HPCExecutor(
            username=username,
            address=address,
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

        1. The `workflow_executor` keyword argument can be removed once [Issue 1024](https://github.com/Quantum-Accelerators/quacc/issues/1024) is resolved.

        !!! Hint

            The most common cause of issues is related to the job scheduler details (i.e. the `resource_spec_kwargs` and the `job_attributes_kwargs`). If your job fails on the remote machine, check the `~/.psij` directory for a history and various log files associated with your attempted job submissions.

    === "Tiger"

        Coming soon.

=== "Parsl ⭐"

    From the login node of the remote machine, run the following once:

    === "Perlmutter"

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

    === "Tiger"

        Coming soon.

    Then run the following on the login node of the remote machine:

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

## Example 2: ORCA

=== "Covalent ⭐"

=== "Parsl ⭐"

    Coming soon.

=== "Jobflow"

    Coming soon.

## Example 3: VASP

=== "Covalent ⭐"

=== "Parsl ⭐"

    Coming soon.

=== "Jobflow"

    Coming soon.
