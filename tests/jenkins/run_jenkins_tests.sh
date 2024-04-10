#!/bin/bash

# ------- Path Setup -------
source ~/.bashrc

# ------- Package Setup -------
if ! conda env list | grep 'quacc' > /dev/null; then
    conda create --name quacc python=3.11
conda activate quacc
pip install uv
uv pip install --force-reinstall --no-deps "ase @ https://gitlab.com/ase/ase/-/archive/master/ase-master.zip"
uv pip install -r tests/requirements.txt "quacc[dev] @ ."

# ------- Request and run Slurm job -------
salloc -N 1 -n 32 -t 00:10:00 tests/jenkins/hpc_tests.sh
