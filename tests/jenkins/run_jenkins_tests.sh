#!/bin/bash

# Setup
source ~/.bashrc
conda activate quacc
pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
pip install -e .[dev]

# Gaussian, ORCA, GULP, VASP
module load gaussian/g16
module load openmpi/gcc/4.1.2
salloc -N 1 -n 32 -t 00:10:00 pytest \
    tests/local/recipes/gaussian_recipes \
    tests/local/recipes/gulp_recipes \
    tests/local/recipes/orca_recipes \
    tests/local/recipes/vasp_recipes/jenkins \
    --noconftest
