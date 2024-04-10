#!/bin/bash

# Setup
source ~/.bashrc
conda activate quacc
pip install uv
uv pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
uv pip install -r tests/requirements.txt "quacc[dev] @ ."

# Request and run Slurm job
salloc -N 1 -n 32 -t 00:10:00 tests/jenkins/hpc_tests.sh
