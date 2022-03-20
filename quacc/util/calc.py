import os
from copy import deepcopy

from ase.atoms import Atoms
from monty.tempfile import ScratchDir

from quacc import SETTINGS


def run_calc(
    atoms: Atoms,
    scratch_dir: str = SETTINGS.SCRATCH_DIR,
    gzip: bool = SETTINGS.GZIP_FILES,
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
    scratch_dir : str
        Path where a tmpdir should be made for running the calculation. If None,
        the current working directory will be used.
    gzip : bool
        Whether to gzip the output files.
    copy_from_store_dir : bool
        Whether to copy any pre-existing files from the original directory to the
        scratch_dir before running the calculation.

    Returns
    -------
    .Atoms
        The updated .Atoms object,
    """

    atoms = deepcopy(atoms)
    scratch_dir = scratch_dir or os.getcwd()
    if atoms.calc is None:
        raise ValueError("Atoms object must have attached calculator.")

    with ScratchDir(
        os.path.abspath(scratch_dir),
        create_symbolic_link=os.name != "nt",
        copy_from_current_on_enter=copy_from_store_dir,
        copy_to_current_on_exit=True,
        gzip_on_exit=gzip,
        delete_removed_files=False,
    ):

        # Run calculation via get_potential_energy()
        atoms.get_potential_energy()

    return atoms
