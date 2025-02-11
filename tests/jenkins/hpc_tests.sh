#!/bin/bash

# Fail qucickly
set -e

# Define path to directories on Adroit
export SCRATCH=/scratch/network/ROSENGROUP
export SOFTWARE=$SCRATCH/software

# Ensure results are written to scratch directory
export QUACC_RESULTS_DIR=$SCRATCH/jenkins

# Set env vars
export GULP_LIB=$SOFTWARE/gulp/gulp-6.1.2/Libraries
export PATH=$SOFTWARE/gulp/gulp-6.1.2/bin:$SOFTWARE/orca_6_0_0_shared_openmpi416:$SOFTWARE/vasp/vasp.6.5.0/bin:$PATH
export QUACC_VASP_PARALLEL_CMD="srun -N 1 --ntasks-per-node 32"
export VASP_PP_PATH=$SOFTWARE/vasp
export ASE_VASP_VDW=$SOFTWARE/vasp/vdw_kernel

# GULP
module purge
pytest tests/core/recipes/gulp_recipes --noconftest

# Gaussian
module purge
module load gaussian/g16
pytest tests/core/recipes/gaussian_recipes --noconftest

# ORCA
module purge
module load nvhpc/24.5 openmpi/nvhpc-24.5/4.1.6
pytest tests/core/recipes/orca_recipes --noconftest

# VASP
module purge
module load intel/2021.1.2 intel-mpi/intel/2021.3.1 intel-mkl/2021.1.1
pytest tests/core/recipes/vasp_recipes/jenkins --noconftest

# Q-Chem
module purge
pytest tests/core/recipes/qchem_recipes/jenkins --noconftest

# Clean up
rm -r $QUACC_RESULTS_DIR
