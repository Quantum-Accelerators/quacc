# Installation

## Pre-Requisites

To use QuAcc, you must have a Python environment set up. For this, we recommend using an environment manager, such as Miniconda. To install Python via Miniconda, download and run your operating system's [Miniconda installer](https://docs.conda.io/en/latest/miniconda.html).

With Miniconda installed, create a fresh environment by running `conda create -n quacc python=3.9`. You can then activate this environment via `conda activate quacc`.

## Installing QuAcc

In your newly activated conda environment, you can install QuAcc and its related dependencies by using `pip install quacc`. For the development version of QuAcc (currently recommended), you can instead use `pip install git+https://github.com/arosen93/quacc.git`.

QuAcc can be installed with several `extras`, as outlined in the `setup.py` file. These are described below:
- `quacc[fireworks]`: Installs dependencies to enable the use of Jobflow and FireWorks.
- `quacc[vasp]`: Installs dependencies to enable the use of VASP with Custodian.
- `quacc[xtb]`: Installs dependencies to enable the use of xtb via tblite.
- `quacc[dev]`: Installs dependencies to enable local development of QuAcc.

If you plan to use QuAcc with a workflow manager, ensure that you install QuAcc on all machines where the code will run.