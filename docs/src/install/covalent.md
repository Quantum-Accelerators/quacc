# Covalent Setup

## Introduction

QuAcc is made to directly interface with [Covalent](https://github.com/AgnostiqHQ/covalent), an open-source workflow management and dispatching solution. In general, we recommend using Covalent to run QuAcc calculations, although it is not strictly required in most cases.

We will describe only the most important setup instructions here. For full details, refer to the [Covalent documentation](https://covalent.readthedocs.io/en/latest/index.html). Please refer to their [Discussion Board](https://github.com/AgnostiqHQ/covalent/discussions) for any Covalent-specific questions.

## Installation

Covalent was installed when you installed QuAcc, so that part is already taken care of. Depending on where you wish to run your QuAcc calculations, you may also need to install the relevant plugin for the type of machine you wish to execture your calculations on, as described [here](https://covalent.readthedocs.io/en/latest/plugins.html).

Note that you don't need to install a specific plugin to run QuAcc locally. So, for testing purposes, you can skip this section. However, for production-quality calculations, we anticipate that most users will use Slurm. To do so, you will need to install the [Slurm executor plugin](https://covalent.readthedocs.io/en/latest/api/executors/slurm.html) via `pip install covalent-slurm-plugin`. Refer to the docuementation for the Slurm executor plugin for how to configure it for your specific high-performance computing machine (or machines!) of choice. 

## Setting the Covalent Configuration

There are two Covalent configuration settings you must modify from their default values to get QuAcc to work properly. Read [here](https://covalent.readthedocs.io/en/latest/how_to/config/customization.html) for how to customize your configuration settings.

First, you must ensure that Covalent runs every calculation in a unique working directory otherwise parallel jobs are likely to crash, as most quantum chemistry codes write out files to the current working directory by default. To address this, set `create_unique_workdir=True` within the options of each executor that you plan to use.

Second, if you are using the Slurm executor, you must set `use_srun=False` in your Covalent configuration file (under `executors.slurm.options`) otherwise there will be parallelization issues at runtime.

Here is how we recommend addressing this:

1. Open `~/.config/covalent/covalent.conf`.
2. Under the `[executors.local]` and `[executors.dask]` sections, set `create_unique_workdir = True`.
3. Under the `[executors.slurm.options]` section (if present), set `use_srun = False`. If this section is not present but you plan to use Slurm, it will appear after you [install the Slurm plugin](https://covalent.readthedocs.io/en/latest/api/executors/slurm.html#installation).
4. Save your edits to this file.

## Starting the Server

Covalent uses a server to store calculation details and results. To start the server, simply run `covalent start` in your terminal when you wish to run calculations. Once you start scaling up your calculations, we recommend hosting the Covalent server on a dedicated machine, but for now you can simply run it on your local machine.