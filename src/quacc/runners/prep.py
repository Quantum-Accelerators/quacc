"""Prepration for runners."""
from __future__ import annotations

import os
from pathlib import Path
from shutil import move, rmtree
from typing import TYPE_CHECKING

from monty.shutil import gzip_dir, remove

from quacc import SETTINGS
from quacc.utils.files import (
    copy_decompress_files,
    copy_decompress_tree,
    make_unique_dir,
)

if TYPE_CHECKING:
    from ase.atoms import Atoms


def calc_setup(
    atoms: Atoms, copy_files: str | Path | list[str | Path] | None = None
) -> tuple[Path, Path]:
    """
    Perform staging operations for a calculation, including copying files to the scratch
    directory, setting the calculator's directory, decompressing files, and creating a
    symlink to the scratch directory.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on. Must have a calculator
        attached.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    Path
        The path to the unique tmpdir, where the calculation will be run. It will be
        deleted after the calculation is complete. By default, this will be
        located within the `SETTINGS.SCRATCH_DIR`, but if that is not set, it will
        be located within the `SETTINGS.RESULTS_DIR`. For conenience, a symlink
        to this directory will be made in the `SETTINGS.RESULTS_DIR`.
    Path
        The path to the results_dir, where the files will ultimately be stored.
        By defualt, this will be the `SETTINGS.RESULTS_DIR`, but if
        `SETTINGS.CREATE_UNIQUE_DIR` is set, it will be a unique directory
        within the `SETTINGS.RESULTS_DIR`.
    """

    # Create a tmpdir for the calculation
    tmpdir_base = SETTINGS.SCRATCH_DIR or SETTINGS.RESULTS_DIR
    tmpdir = make_unique_dir(base_path=tmpdir_base, prefix="tmp-quacc-")

    # Set the calculator's directory
    atoms.calc.directory = tmpdir

    # Define the results directory
    job_results_dir = SETTINGS.RESULTS_DIR
    if SETTINGS.CREATE_UNIQUE_DIR:
        job_results_dir /= f"{tmpdir.name.split('tmp-')[-1]}"

    # Create a symlink to the tmpdir
    if os.name != "nt" and SETTINGS.SCRATCH_DIR:
        symlink = SETTINGS.RESULTS_DIR / f"symlink-{tmpdir.name}"
        symlink.unlink(missing_ok=True)
        symlink.symlink_to(tmpdir, target_is_directory=True)

    # Copy files to tmpdir and decompress them if needed
    if isinstance(copy_files, dict):
        copy_decompress_tree(copy_files, tmpdir)
    if isinstance(copy_files, list):
        copy_decompress_files(copy_files, tmpdir)
    elif isinstance(copy_files, (str, Path)):
        copy_decompress_files([copy_files], tmpdir)

    # NOTE: Technically, this breaks thread-safety since it will change the cwd
    # for all threads in the current process. However, elsewhere in the code,
    # we use absolute paths to avoid issues. We keep this here for now because some
    # old ASE calculators do not support the `directory` keyword argument.
    os.chdir(tmpdir)

    return tmpdir, job_results_dir


def calc_cleanup(atoms: Atoms, tmpdir: Path | str, job_results_dir: Path | str) -> None:
    """
    Perform cleanup operations for a calculation, including gzipping files, copying
    files back to the original directory, and removing the tmpdir.

    Parameters
    ----------
    atoms
        The Atoms object after the calculation. Must have a calculator
        attached.
    tmpdir
        The path to the tmpdir, where the calculation will be run. It will be
        deleted after the calculation is complete.
    job_results_dir
        The path to the job_results_dir, where the files will ultimately be
        stored. A symlink to the tmpdir will be made here during the calculation
        for convenience.

    Returns
    -------
    None
    """

    job_results_dir, tmpdir = Path(job_results_dir), Path(tmpdir)

    # Reset the calculator's directory
    atoms.calc.directory = job_results_dir

    # Make the results directory
    job_results_dir.mkdir(parents=True, exist_ok=True)

    # NOTE: Technically, this breaks thread-safety since it will change the cwd
    # for all threads in the current process. However, elsewhere in the code,
    # we use absolute paths to avoid issues. We keep this here for now because some
    # old ASE calculators do not support the `directory` keyword argument.
    os.chdir(job_results_dir)

    # Gzip files in tmpdir
    if SETTINGS.GZIP_FILES:
        gzip_dir(tmpdir)

    # Move files from tmpdir to job_results_dir
    for file_name in os.listdir(tmpdir):
        if Path(job_results_dir / file_name).is_dir():
            remove(job_results_dir / file_name)
        move(tmpdir / file_name, job_results_dir / file_name)

    # Remove symlink to tmpdir
    if os.name != "nt" and SETTINGS.SCRATCH_DIR:
        symlink_path = SETTINGS.RESULTS_DIR / f"symlink-{tmpdir.name}"
        symlink_path.unlink(missing_ok=True)

    # Remove the tmpdir
    rmtree(tmpdir, ignore_errors=True)
