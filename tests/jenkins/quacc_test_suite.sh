source ~/.bashrc
conda activate quacc
pip install --force-reinstall --no-deps https://gitlab.com/ase/ase/-/archive/master/ase-master.zip
pip install -e .[dev]

module load gaussian/g16 # gaussian
module load openmpi/gcc/4.1.2 # orca
salloc -N 1 -n 4 -t 00:10:00 pytest tests/local/recipes/gaussian_recipes tests/local/recipes/gulp_recipes tests/local/recipes/orca_recipes --noconftest --cov=quacc --cov-report=xml

# module purge
# module load intel/2021.1.2 # vasp
# module load intel-mpi/intel/2021.3.1 # vasp
# module load hdf5/intel-2021.1/1.10.6 # vasp
# salloc -N 1 -n 32 -t 00:10:00 pytest tests/local/recipes/vasp_recipes --noconftest --cov=quacc --cov-report=xml

codecovcli upload-process
