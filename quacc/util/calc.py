import os
from copy import deepcopy
from shutil import copy, rmtree
from tempfile import mkdtemp

from ase.atoms import Atoms
from monty.shutil import copy_r, gzip_dir


def run_calc(
    atoms: Atoms, store_dir: str = None, scratch_basedir: str = None, gzip: bool = False
) -> float:
    """
    Run a calculation in a scratch directory and copy the results back to the
    original directory. This can be useful if file I/O is slow in the working
    directory, so long as file transfer speeds are reasonable.

    This is a wrapper around atoms.get_potential_energy().

    Parameters
    ----------
    atoms : .Atoms
        The Atoms object to run the calculation on.
    store_dir : str
        Path where calculation results should be stored. Also will copy all files
        from this directory at runtime. If None, the current working directory will be used.
    scratch_basedir : str
        Path to the base directory where a tmp directory will be made for
        scratch files. If None, a temporary directory in $SCRATCH will be used.
        If $SCRATCH is not present, a tmp directory will be made in store_dir.
    gzip : bool
        Whether to gzip the output files.

    Returns
    -------
    .Atoms
        The updated .Atoms object,
    """

    atoms = deepcopy(atoms)
    if atoms.calc is None:
        raise ValueError("Atoms object must have attached calculator.")

    # Find the relevant paths
    if not store_dir:
        store_dir = os.getcwd()
    if not scratch_basedir:
        if "SCRATCH" in os.environ:
            scratch_basedir = os.path.expandvars("$SCRATCH")
        else:
            scratch_basedir = store_dir
    if not os.path.exists(scratch_basedir):
        raise OSError(f"Cannot find {scratch_basedir}")

    scratch_path = mkdtemp(dir=scratch_basedir, prefix="quacc-")

    # Copy files to scratch
    for f in os.listdir(store_dir):
        if os.path.isfile(os.path.join(store_dir, f)):
            copy(os.path.join(store_dir, f), os.path.join(scratch_path, f))

    # Leave a note in the run directory for where the scratch is located in case
    # the job dies partway through
    scratch_path_note = os.path.join(store_dir, "scratch_path.txt")
    with open(scratch_path_note, "w") as f:
        f.write(scratch_path)

    # Run calculation via get_potential_energy()
    atoms.calc.directory = scratch_path
    atoms.get_potential_energy()

    # gzip files and recursively copy files back to store_dir
    if gzip:
        gzip_dir(scratch_path)
    copy_r(scratch_path, store_dir)

    # Remove the scratch note and scratch dir
    if os.path.exists(scratch_path_note):
        os.remove(scratch_path_note)
    if os.path.exists(scratch_path):
        rmtree(scratch_path)

    return atoms
