#!/bin/bash

# Setup
source ~/.bashrc
conda activate quacc
pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
pip install -e .[dev]

# Gaussian, ORCA, GULP
module load gaussian/g16
module load openmpi/gcc/4.1.2
salloc -N 1 -n 32 -t 00:10:00 pytest \
    tests/core/recipes/gaussian_recipes \
    tests/core/recipes/gulp_recipes \
    tests/core/recipes/orca_recipes \
    --noconftest

# VASP
module purge
module load intel/2021.1.2
module load intel-mpi/intel/2021.3.1
module load hdf5/intel-2021.1/1.10.6
salloc -N 1 -n 32 -t 00:10:00 pytest tests/core/recipes/vasp_recipes/jenkins --noconftest

# Upload to codecov
codecovcli upload-process
