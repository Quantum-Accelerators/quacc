import os
from shutil import copy, copyfileobj
from ase.atoms import Atoms


def run_calc(
    atoms: Atoms, run_path: str = None, scratch_path: str = None, gzip: bool = False
) -> None:
    """
    Run a calculation in a scratch directory and copy the results back to the
    original directory. This can be useful if file I/O is slow in the working
    directory, so long as file transfer speeds are reasonable.

    This is a wrapper around atoms.get_potential_energy().

    Parameters
    ----------
    atoms : .Atoms
        The Atoms object to run the calculation on.
    run_path : str
        Path to the directory containing the calculation to be run.
        If None, the current working directory will be used.
    scratch_path : str
        Path to the directory to be used as the scratch directory.
        If None, a temporary directory in $SCRATCH will be used.
    gzip : bool
        Whether to gzip the output files.
    """

    # Find the relevant paths
    if not run_path:
        run_path = os.getcwd()
    if not scratch_path:
        scratch_path = os.path.expandvars("$SCRATCH")
        if not scratch_path:
            raise EnvironmentError(
                "scratch_path is None yet $SCRATCH environment variable is not set"
            )

    # Copy files from working directory to scratch directory
    for f in os.listdir(run_path):
        copy(os.path.join(run_path, f), os.path.join(scratch_path, f))

    # Run calculation via get_potential_energy()
    os.chdir(scratch_path)
    atoms.get_potential_energy()
    os.chdir(run_path)

    # Copy files from scratch directory to working directory
    for f in os.listdir(scratch_path):
        # gzip the files in scratch before copying them back
        if gzip:
            with open(os.path.join(scratch_path, f), "rb") as f_in:
                with gzip.open(os.path.join(run_path, f + ".gz"), "wb") as f_out:
                    copyfileobj(f_in, f_out)
        else:
            copy(os.path.join(scratch_path, f), os.path.join(run_path, f))
