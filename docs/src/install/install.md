# Installation

## Pre-Requisites

```{note}
If you are on a Windows machine, you will need to install and use the [Windows Subsystem for Linux (WSL)](https://ubuntu.com/wsl).
```

1. Install Python, preferably in conjunction with an environment manager. For instance, download and run the [Miniconda installer](https://docs.conda.io/en/latest/miniconda.html).
2. When asked to add Miniconda to your `PATH`, select yes.
3. Create a fresh Python environment by running `conda create -n quacc python=3.9`.
4. Activate this environment via `conda activate quacc`.

## Installing Quacc

```{note}
If you plan to use Quacc with a workflow manager, ensure that you install Quacc on all machines where the code will run.
```

In your newly activated conda environment, you can install Quacc and its related dependencies by running `pip install quacc`. For the development version of Quacc (**currently required**), you should instead use `pip install git+https://github.com/arosen93/quacc.git`.

Quacc can be installed with several "extras," as outlined in the `setup.py` file. These are described below:

- `quacc[fireworks]`: Installs dependencies to enable the use of FireWorks.
- `quacc[tblite]`: Installs dependencies to enable the use of tblite.
- `quacc[dev]`: Installs dependencies to enable local development of Quacc.
