#!/bin/bash

# ------- Path Setup -------
source ~/.bashrc

# ------- Package Setup -------
module load anaconda3/2024.10
conda activate quacc
conda install -c conda-forge openbabel
pip install uv
uv pip install -r tests/requirements.txt "quacc[dev] @ ."

# ------- Request and run Slurm job -------
salloc -N 1 -n 32 -t 00:15:00 tests/jenkins/hpc_tests.sh
