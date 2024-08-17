#!/bin/bash

# Fail fast
set -e

# Define path to directories on Adroit
export SCRATCH=/scratch/network/ROSENGROUP
export SOFTWARE=$SCRATCH/software

# Ensure results are written to scratch directory
export QUACC_RESULTS_DIR=$SCRATCH/jenkins

# GULP
export GULP_LIB=$SOFTWARE/gulp/gulp-6.1.2/Libraries
export PATH=$SOFTWARE/gulp/gulp-6.1.2/bin:$PATH
module purge
pytest tests/core/recipes/gulp_recipes --noconftest

# Gaussian
module purge
module load gaussian/g16
pytest tests/core/recipes/gaussian_recipes --noconftest

# ORCA
export PATH=$SOFTWARE/orca:$PATH
module purge
module load openmpi/gcc/4.1.2
pytest tests/core/recipes/orca_recipes --noconftest

# VASP
export QUACC_VASP_PARALLEL_CMD="srun -N 1 --ntasks-per-node 32"
export VASP_PP_PATH=$SOFTWARE/vasp
export ASE_VASP_VDW=$SOFTWARE/vasp/vdw_kernel.bindat
export PATH=$SOFTWARE/vasp/vasp.6.4.3/bin:$PATH
module purge
module load intel/2021.1.2 intel-mpi/intel/2021.3.1 hdf5/intel-2021.1/1.10.6
pytest tests/core/recipes/vasp_recipes/jenkins --noconftest

# Q-Chem
module purge
pytest tests/core/recipes/qchem_recipes/jenkins --noconftest

# Clean up
rm -r $QUACC_RESULTS_DIR
