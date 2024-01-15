#!/bin/bash

# Setup
source ~/.bashrc
conda activate quacc
pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
pip install -e .[dev]

# Gaussian, ORCA, GULP
salloc -N 1 -n 32 -t 00:10:00 \
    module purge && pytest tests/core/recipes/gulp_recipes --noconftest \
    module purge && module load gaussian/g16 && pytest tests/core/recipes/gaussian_recipes --noconftest \
    module purge && module load openmpi/gcc/4.1.2 && pytest tests/core/recipes/orca_recipes --noconftest \
    module purge && module load intel/2021.1.2 intel-mpi/intel/2021.3.1 hdf5/intel-2021.1/1.10.6 && pytest tests/core/recipes/vasp_recipes/jenkins --noconftest
