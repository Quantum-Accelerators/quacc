# Installation

## Pre-Requisites

!!! note

    If you are on a Windows machine, you will need to install and use the [Windows Subsystem for Linux (WSL)](https://ubuntu.com/wsl).

1. Install Python, preferably in conjunction with an environment manager. For instance, download and run the [Miniconda installer](https://docs.conda.io/en/latest/miniconda.html)
2. When asked to add Miniconda to your `PATH`, select yes
3. Create a fresh Python environment by running `conda create -n quacc python=3.10`
4. Activate this environment via `conda activate quacc`

## Installing quacc

In your newly activated conda environment, run the following commands to install quacc. Note that you will need to install quacc on all machines where you plan to run calculations.

```bash
# Install development version of quacc
pip install git+https://github.com/quantum-accelerators/quacc.git

# Set default configuration parameters
quacc config
```

!!! Tip

    Everything beyond this point in the installation guide is to add on useful features to quacc. So, if you are just getting started, check out the [Quacc Basics](../user/basics.md) page. Then come back to installing additional features as you need them.

## Optional Dependencies

Quacc can be installed with several "extras," as outlined in the `pyproject.toml` file. These are described below.

### Calculators

- `quacc[tblite]`: Installs dependencies to enable the use of tblite.

### Workflow Managers

- `quacc[fireworks]`: Installs dependencies to enable the use of FireWorks.
- `quacc[parsl]`: Installs dependencies to enable the use of Parsl.

### Miscellaneous

- `quacc[defects]`: Installs dependencies to enable the use of defect workflows.
- `quacc[optimizers]`: Installs dependencies to enable the use of the Sella optimizer.

### Development

- `quacc[dev]`: Installs dependencies to enable local development of quacc.
- `quacc[docs]`: Installs dependencies to build the documentation.
- `quacc[strict]`: Installs dependencies that match the test suite on GitHub.
