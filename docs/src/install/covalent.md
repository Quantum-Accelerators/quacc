# Covalent Setup

## Introduction

Quacc is made to directly interface with [Covalent](https://github.com/AgnostiqHQ/covalent), an open-source workflow management and dispatching solution. We will describe only the most important setup instructions here. For full details, refer to the [Covalent documentation](https://covalent.readthedocs.io/en/latest/index.html). Please refer to the Covalent [Discussion Board](https://github.com/AgnostiqHQ/covalent/discussions) for any Covalent-specific questions.


```{note}
If you plan to use Jobflow as your workflow manager instead of Covalent, you may skip this section.
```

## Starting the Server

Covalent uses a server to dispatch and retrieve calculation details and results. To start the server, simply run `covalent start` in your terminal when you wish to run calculations.

Once you start scaling up your calculations, we recommend hosting the Covalent server on a dedicated machine or using [Covalent Cloud](https://www.covalent.xyz/cloud/), but for now you can simply run it on your local machine.

## Plugin Installation

```{note}
If you are just starting out, try running some test calculations locally first. Then come back and install the relevant Covalent plugins for your machines.
```

Depending on where you wish to run your Quacc calculations, you may need to install the relevant Covalent plugin for the type of machine you wish to execute your calculations on, as described in the [Covalent plugin documention](https://covalent.readthedocs.io/en/latest/plugins.html).

Note that you don't need to install a specific plugin to run Quacc locally. So, if you are new to Quacc, you should skip this section and come back later. However, for production-quality calculations, we anticipate that most users will use Slurm. To do so, you will need to install the [Slurm executor plugin](https://covalent.readthedocs.io/en/latest/api/executors/slurm.html) via `pip install covalent-slurm-plugin`. Refer to the docuementation for the Slurm executor plugin for how to configure it for your specific high-performance computing machine (or machines!) of choice.
