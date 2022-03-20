import os
from copy import deepcopy

from ase.atoms import Atoms
from monty.tempfile import ScratchDir

from quacc import SETTINGS


def run_calc(
    atoms: Atoms,
    scratch_dir: str = SETTINGS.SCRATCH_DIR,
    gzip: bool = True,
    copy_from_store_dir: bool = False,
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
        Path where calculation results should be stored.
        If None, the current working directory will be used.
    scratch_dir : str
        Path to the base directory where a tmp directory will be made for
        scratch files. If None, a temporary directory in $SCRATCH will be used.
        If $SCRATCH is not present, the current working directory will be used.
    gzip : bool
        Whether to gzip the output files.
    copy_from_store_dir : bool
        Whether to copy any pre-existing files from the store_dir to the scratch_dir
        before running the calculation.

    Returns
    -------
    .Atoms
        The updated .Atoms object,
    """

    atoms = deepcopy(atoms)
    if atoms.calc is None:
        raise ValueError("Atoms object must have attached calculator.")

    # Find the relevant paths
    if not scratch_dir:
        if "SCRATCH" in os.environ:
            scratch_dir = os.path.expandvars("$SCRATCH")
        else:
            scratch_dir = os.getcwd()

    if not os.path.exists(scratch_dir):
        raise ValueError(f"Cannot find {scratch_dir}")

    with ScratchDir(
        os.path.abspath(scratch_dir),
        create_symbolic_link=False if os.name == "nt" else True,
        copy_from_current_on_enter=copy_from_store_dir,
        copy_to_current_on_exit=True,
        gzip_on_exit=gzip,
        delete_removed_files=False,
    ):

        # Run calculation via get_potential_energy()
        atoms.get_potential_energy()

    return atoms
