import os
from copy import deepcopy
import gzip
from shutil import copy, copytree, copyfileobj
from tempfile import TemporaryDirectory
from ase.atoms import Atoms


def run_calc(
    atoms: Atoms, run_dir: str = None, scratch_dir: str = None, compress: bool = False
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
    run_dir : str
        Path to the directory containing the calculation to be run.
        If None, the current working directory will be used.
    scratch_dir : str
        Path to the base directory to store the scratch temp directories.
        If None, a temporary directory in $SCRATCH will be used. If $SCRATCH
        is not present, everything will be run in run_dir.
    compress : bool
        Whether to compress the output files.

    Returns
    -------
    .Atoms
        The updated .Atoms object,
    """

    atoms = deepcopy(atoms)
    if atoms.calc is None:
        raise ValueError("Atoms object must have attached calculator.")

    # Find the relevant paths
    if not run_dir:
        run_dir = os.getcwd()
    if not scratch_dir:
        if "SCRATCH" in os.environ:
            scratch_dir = os.path.expandvars("$SCRATCH")
        else:
            scratch_dir = run_dir

    with TemporaryDirectory(dir=scratch_dir, prefix="quacc-") as scratch_path:

        _copy_to_scratch(run_dir, scratch_path)

        # Leave a note in the run directory for where the scratch is located in case
        # the job dies partway through
        scratch_path_note = os.path.join(run_dir, "scratch_path.txt")
        with open(scratch_path_note, "w") as f:
            f.write(scratch_path)

        # Run calculation via get_potential_energy()
        os.chdir(scratch_path)
        atoms.get_potential_energy()
        os.chdir(run_dir)

        _copy_from_scratch(run_dir, scratch_path, compress=compress)

    # Remove the scratch note
    if os.path.exists(scratch_path_note):
        os.remove(scratch_path_note)

    return atoms


def _copy_to_scratch(run_dir: str, scratch_path: str):
    """
    Copy files from working directory to scratch directory.

    Parameters
    ----------
    run_dir : str
        Path to the run directory
    scratch_path : str
        Path to the scratch directory
    """
    # Copy files from working directory to scratch directory
    for f in os.listdir(run_dir):
        if os.path.isdir(os.path.join(run_dir, f)):
            copytree(os.path.join(run_dir, f), os.path.join(scratch_path, f))
        else:
            copy(os.path.join(run_dir, f), os.path.join(scratch_path, f))


def _copy_from_scratch(run_dir: str, scratch_path: str, compress: bool = False):
    """
    Copy files from scratch directory to working directory.

    Parameters
    ----------
    run_dir : str
        Path to the run directory
    scratch_path : str
        Path to the scratch directory
    """
    # Copy files from scratch directory to working directory
    for f in os.listdir(scratch_path):
        if os.path.isdir(os.path.join(scratch_path, f)):
            copytree(os.path.join(scratch_path, f), os.path.join(run_dir, f))
            # NOTE: Would be worth .tar.gz'ing the directory if compress is True.
        else:
            # gzip the files in scratch before copying them back
            if compress and os.path.splitext(f)[1] != ".gz":
                with open(os.path.join(scratch_path, f), "rb") as f_in:
                    with gzip.open(os.path.join(run_dir, f + ".gz"), "wb") as f_out:
                        copyfileobj(f_in, f_out)
            else:
                copy(os.path.join(scratch_path, f), os.path.join(run_dir, f))
