# Covalent Setup

## Introduction

Quacc is made to directly interface with [Covalent](https://github.com/AgnostiqHQ/covalent), an open-source workflow management and dispatching solution with a convenient user interface for job monitoring.

```{note}
If you plan to use Parsl or Jobflow as your workflow manager instead of Covalent, skip to the [Parsl Setup](advanced/parsl.md) or [Jobflow Setup](advanced/jobflow.md) sections, respectively.
```

## Starting the Server

Covalent uses a server to dispatch and retrieve calculation details and results. To start the server, simply run `covalent start` in your terminal. It will return a URL (typically http://localhost:48008) that you can use to access the Covalent dashboard, which is shown below.

![Covalent UI](../_static/install/ui_blank.jpg)

```{tip}
Once you start scaling up your calculations, we recommend hosting the Covalent server on a dedicated machine or in the cloud. Refer to the [Covalent Deployment Guide](https://docs.covalent.xyz/docs/user-documentation/server-deployment) for details.
```

## Plugin Installation

```{note}
If you are just starting out, try running some test calculations locally first. Then come back and install the relevant Covalent plugins for your machines.
```

Depending on where you wish to run your Quacc calculations, you may need to install the corresponding Covalent plugin, as described in the [Covalent plugin documentation](https://docs.covalent.xyz/docs/features/executor-plugins/exe). For production-quality calculations, we anticipate that most users will rely on the `SlurmExecutor`, which can be installed via `pip install covalent-slurm-plugin`. Refer to the [Slurm executor documentation](https://docs.covalent.xyz/docs/user-documentation/api-reference/executors/slurm) for instructions on how to configure it for your desired high-performance computing machine.

For submitting jobs to Perlmutter at NERSC from your local machine, an example `SlurmExecutor` configuration with support for an [`sshproxy`](https://docs.nersc.gov/connect/mfa/#sshproxy)-based multi-factor authentication certificate might look like the following:

```python
executor = ct.executor.SlurmExecutor(
    username="YourUserName",
    address="perlmutter-p1.nersc.gov",
    ssh_key_file="~/.ssh/nersc",
    cert_file="~/.ssh/nersc-cert.pub",
    remote_workdir="$SCRATCH",
    conda_env="quacc",
    options={
        "nodes": 1,
        "qos": "debug",
        "constraint": "cpu",
        "account": "YourAccountName",
        "job-name": "quacc",
        "time": "00:10:00",
    },
    prerun_commands=[
        "export COVALENT_CONFIG_DIR=$SCRATCH",
        "export OMP_PROC_BIND=spread",
        "export OMP_PLACES=threads",
        "export OMP_NUM_THREADS=1",
    ],
)
```
