# Installation

## Pre-Requisites

To use QuAcc, you must have a Python environment set up. For this, we recommend using an environment manager, such as Miniconda. To install Python via Miniconda, download and run your operating system's [Miniconda installer](https://docs.conda.io/en/latest/miniconda.html).

With Miniconda installed, create a fresh environment by running `conda create -n quacc python=3.9`. You can then activate this environment via `conda activate quacc`.

## Installing QuAcc

In your new activated conda environment, you can install QuAcc and its related dependencies by using `pip install quacc`. For the development version of QuAcc (recommended), you can instead use `pip install git+https://github.com/arosen93/quacc.git`.

To install QuAcc with support for all features out-of-the-box, append `[all]` to the end (e.g. `pip install git+https://github.com/arosen93/quacc.git`).
