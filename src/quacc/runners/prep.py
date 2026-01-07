"""Prepration for runners."""

from __future__ import annotations

import os
from logging import getLogger
from pathlib import Path
from shutil import move, rmtree
from typing import TYPE_CHECKING

from monty.shutil import gzip_dir

from quacc import JobFailure, get_settings
from quacc.utils.files import copy_decompress_files, make_unique_dir

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import Filenames, SourceDirectory

LOGGER = getLogger(__name__)


def calc_setup(
    atoms: Atoms | None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> tuple[Path, Path]:
    """
    Perform staging operations for a calculation, including copying files to the scratch
    directory, setting the calculator's directory, decompressing files, and creating a
    symlink to the scratch directory.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on. Must have a calculator
        attached. If None, no modifications to the calculator's directory will be made.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    Path
        The path to the unique tmpdir, where the calculation will be run. It will be
        deleted after the calculation is complete. By default, this will be
        located within `QuaccSettings.SCRATCH_DIR`, but if that is not set, it will
        be located within the `QuaccSettings.RESULTS_DIR`. For conenience, a symlink
        to this directory will be made in the `QuaccSettings.RESULTS_DIR`.
    Path
        The path to the results_dir, where the files will ultimately be stored.
        By defualt, this will be the `QuaccSettings.RESULTS_DIR`, but if
        `QuaccSettings.CREATE_UNIQUE_DIR` is set, it will be a unique directory
        within the `QuaccSettings.RESULTS_DIR`.
    """
    settings = get_settings()
    base_path = settings.SCRATCH_DIR or settings.RESULTS_DIR

    # If CREATE_UNIQUE_DIR is True, make a unique tmpdir for the calculation
    if settings.CREATE_UNIQUE_DIR:
        rundir = make_unique_dir(base_path=base_path, prefix="tmp-quacc-")
        job_results_dir = settings.RESULTS_DIR
        job_results_dir /= f"{rundir.name.split('tmp-')[-1]}"

        # Create a symlink to the rundir
        if os.name != "nt" and settings.SCRATCH_DIR:
            symlink_path = settings.RESULTS_DIR / f"symlink-{rundir.name}"
            symlink_path.symlink_to(rundir, target_is_directory=True)
    # Else, run in the current working directory
    else:
        rundir = base_path
        job_results_dir = rundir

    LOGGER.info(f"Calculation will run at {rundir}")

    # Set the calculator's directory
    if atoms is not None:
        atoms.calc.directory = rundir

    # Copy files to rundir and decompress them if needed
    if copy_files:
        if isinstance(copy_files, str | Path):
            copy_files = {copy_files: "*"}

        for source_directory, filenames in copy_files.items():
            copy_decompress_files(source_directory, filenames, rundir)

    return rundir, job_results_dir


def calc_cleanup(
    atoms: Atoms | None, rundir: Path | str, job_results_dir: Path | str
) -> None:
    """
    Perform cleanup operations for a calculation, including gzipping files, copying
    files back to the original directory, and removing the tmpdir if made.

    Parameters
    ----------
    atoms
        The Atoms object after the calculation. Must have a calculator
        attached. If None, no modifications to the calculator's directory will be made.
    rundir
        The path to the rundir, where the calculation will be run.
    job_results_dir
        The path to the job_results_dir, where the files will ultimately be
        stored.

    Returns
    -------
    None
    """
    job_results_dir, rundir = Path(job_results_dir), Path(rundir)
    settings = get_settings()

    # Update the calculator's directory
    if atoms is not None:
        atoms.calc.directory = job_results_dir

    # Gzip files in rundir
    if settings.GZIP_FILES:
        gzip_dir(rundir)

    # Safety check before moving files
    if settings.CREATE_UNIQUE_DIR:
        # Safety check
        if "tmp-" not in rundir.name:
            msg = f"{rundir} does not appear to be a tmpdir... exiting for safety!"
            raise ValueError(msg)

        # Move files from rundir to job_results_dir
        move(rundir, job_results_dir)

        # Remove symlink to rundir
        if os.name != "nt" and settings.SCRATCH_DIR:
            symlink_path = settings.RESULTS_DIR / f"symlink-{rundir.name}"
            symlink_path.unlink(missing_ok=True)
    LOGGER.info(f"Calculation results stored at {job_results_dir}")


def terminate(rundir: Path | str, exception: Exception) -> None:
    """
    Terminate a calculation and move files to a failed directory.

    Parameters
    ----------
    rundir
        The path to the rundir, where the calculation was run.
    exception
        The exception that caused the calculation to fail.

    Returns
    -------
    None

    Raises
    -------
    JobFailure
        The exception that caused the calculation to fail plus additional
        metadata.
    """
    rundir = Path(rundir)
    settings = get_settings()

    if settings.CREATE_UNIQUE_DIR:
        job_failed_dir = rundir.with_name(rundir.name.replace("tmp-", "failed-"))
        rundir.rename(job_failed_dir)

        if os.name != "nt" and settings.SCRATCH_DIR:
            old_symlink_path = settings.RESULTS_DIR / f"symlink-{rundir.name}"
            symlink_path = settings.RESULTS_DIR / f"symlink-{job_failed_dir.name}"
            old_symlink_path.unlink(missing_ok=True)
            symlink_path.symlink_to(job_failed_dir, target_is_directory=True)

    else:
        job_failed_dir = rundir

    msg = f"Calculation failed! Files stored at {job_failed_dir}"
    LOGGER.info(msg)

    raise JobFailure(job_failed_dir, message=msg, parent_error=exception) from exception
