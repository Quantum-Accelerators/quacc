# Covalent Setup

## Introduction

Quacc is made to directly interface with [Covalent](https://github.com/AgnostiqHQ/covalent), an open-source workflow management and dispatching solution with a convenient user interface for job monitoring.

```{note}
If you prefer to use a workflow engine other than Covalent, then refer to the ["Alternate Workflow Engines"](alt_workflows/index.md) section of the documentation.
```

## Starting the Server

Covalent uses a server to dispatch and store calculation details and results. To start the server, simply run `covalent start` in your terminal. It will return a URL that you can use to access the Covalent dashboard, which is shown below.

![Covalent UI](../_static/install/ui_blank.jpg)

```{tip}
Once you start scaling up your calculations, we recommend hosting the Covalent server on a dedicated machine or using [Covalent Cloud](https://www.covalent.xyz/cloud/). Refer to the [Covalent Deployment Guide](https://docs.covalent.xyz/docs/user-documentation/server-deployment) for details.
```

## Optional Configuration

Covalent has several [configuration options](https://docs.covalent.xyz/docs/user-documentation/how-to/customization/) that can be modified. Running `quacc config` automatically takes care of setting the ones that are critical for Quacc to run properly.

```{important}
If you ever delete your Covalent configuration (e.g. via `covalent purge`), you will need to re-run `quacc config`.
```

## Plugin Installation

```{note}
If you are just starting out, try running some test calculations locally first. Then come back and install the relevant Covalent plugins for your machines.
```

Depending on where you wish to run your Quacc calculations, you may need to install the corresponding Covalent plugin, as described in the [Covalent plugin documentation](https://docs.covalent.xyz/docs/features/executor-plugins/exe). For production-quality calculations, we anticipate that most users will rely on the `SlurmExecutor`, which can be installed via `pip install covalent-slurm-plugin`.
