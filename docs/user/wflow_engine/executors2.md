# Worked Examples

In this section, we provide a few step-by-step examples of remotely deployed recipes.

!!! Hint

    Before deploying remote calculations for the first time, make sure you can run locally without a workflow manager. Once you are confident that works, try running that same Python script on your desired computing resource (e.g. by submitting it as a job to the scheduler). These preliminary tests will help you identify potential issues early on.

## Pre-Requisites

On the local and remote machines, make a clean Conda environment if you haven't already:

```bash
conda create --name quacc python=3.10
conda activate quacc
```

Then, on both the local and remote machines, install the necessary dependencies:

=== "Covalent ⭐"

    ```bash
    pip install --no-cache-dir https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[covalent]
    quacc set WORKFLOW_ENGINE covalent
    ```

    From the local machine:

    ```bash
    covalent start
    ```

    !!! Note

        If using Perlmutter at NERSC, modify your `~/.bashrc` on the remote machine as follows:

        ```bash title="~/.bashrc"
        export COVALENT_CONFIG_DIR="$SCRATCH/.config/covalent"
        ```

=== "Parsl ⭐"

    ```bash
    pip install --no-cache-dir https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[parsl]
    quacc set WORKFLOW_ENGINE parsl
    ```

=== "Prefect"

    ```bash
    pip install --no-cache-dir https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[prefect]
    quacc set WORKFLOW_ENGINE prefect
    ```

=== "Redun"

    ```bash
    pip install --no-cache-dir https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[redun]
    quacc set WORKFLOW_ENGINE redun
    ```

=== "Jobflow"

    ```bash
    pip install --no-cache-dir https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
    pip install quacc[jobflow]
    quacc set WORKFLOW_ENGINE jobflow
    ```

## Example 1: EMT

When deploying calculations for the first time, it's important to start simple, which is why you should try to run a sample EMT workflow first.

=== "Covalent ⭐"

    ```python
    import covalent as ct
    from ase.build import bulk
    from quacc import flow
    from quacc.recipes.emt.core import relax_job
    from quacc.recipes.emt.slabs import bulk_to_slabs_flow

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
    )

    @flow(executor=executor, workflow_executor=executor) # (1)!
    def workflow(atoms):
        output1 = relax_job(atoms)
        return bulk_to_slabs_flow(output1)

    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id)
    ```

    1. The `workflow_executor` keyword argument can be removed once [Issue 1024](https://github.com/Quantum-Accelerators/quacc/issues/1024) is resolved.

    !!! Hint

        The most common cause of issues is related to the job scheduler details (i.e. the `resource_spec_kwargs` and the `job_attributes_kwargs`). If your job fails on the remote machine, check the `~/.psij` directory for various log files associated with your attempted job submissions. Additionally, set `cleanup=False` in the `HPCExecutor` to retain all PSI/J-related files for debugging purposes.

=== "Parsl ⭐"

    Coming soon.

=== "Prefect"

    Coming soon.

=== "Redun"

    Coming soon.

=== "Jobflow"

    Coming soon.

## Example 2: ORCA

=== "Covalent ⭐"

=== "Parsl ⭐"

    Coming soon.

=== "Prefect"

    Coming soon.

=== "Redun"

    Coming soon.

=== "Jobflow"

    Coming soon.

## Example 4: VASP

=== "Covalent ⭐"

=== "Parsl ⭐"

    Coming soon.

=== "Prefect"

    Coming soon.

=== "Redun"

    Coming soon.

=== "Jobflow"

    Coming soon.
