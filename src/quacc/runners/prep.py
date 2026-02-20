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
from quacc.wflow_tools.context import get_context_path, get_directory_context

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import SourceDirectory
    from quacc.wflow_tools.job_argument import Copy

LOGGER = getLogger(__name__)


def calc_setup(
    atoms: Atoms | None, copy_files: SourceDirectory | Copy | None = None
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
    # Create a tmpdir for the calculation
    settings = get_settings()
    tmpdir_base = (settings.SCRATCH_DIR or settings.RESULTS_DIR).resolve()

    if settings.AUTODISCOVER_DIR:
        # Define the results directory.
        # When AUTODISCOVER_DIR is active, the context module provides both a root
        # directory (directory_context) and a relative path (context_path) that
        # mirrors the flow/subflow/job nesting.
        if Path(get_directory_context()).is_relative_to(settings.RESULTS_DIR):
            # The directory context is already under RESULTS_DIR (normal case).
            job_results_dir = settings.RESULTS_DIR / Path(get_directory_context())
            job_results_dir = job_results_dir / get_context_path()
        else:
            # RESULTS_DIR was overridden (e.g. via swap_settings); use it as the
            # base and append only the context path.
            job_results_dir = settings.RESULTS_DIR / get_context_path()

        # Create a temporary directory with the same internal structure as job results, but with "tmp-" in its
        # top-level directory name. Place it at `tmpdir_base`.
        tmpdir = tmpdir_base / Path(
            "tmp-" + str(job_results_dir.relative_to(settings.RESULTS_DIR))
        )
        tmpdir.mkdir(parents=True, exist_ok=True)
    else:
        tmpdir = make_unique_dir(base_path=tmpdir_base, prefix="tmp-quacc-")
        job_results_dir = settings.RESULTS_DIR.resolve()
        if settings.CREATE_UNIQUE_DIR:
            job_results_dir /= f"{tmpdir.name.split('tmp-')[-1]}"

    LOGGER.info(f"Calculation will run at {tmpdir}")

    # Set the calculator's directory
    if atoms is not None:
        atoms.calc.directory = tmpdir

    # Create a symlink to the tmpdir
    if os.name != "nt" and settings.SCRATCH_DIR:
        symlink_path = settings.RESULTS_DIR / f"symlink-{tmpdir.name}"
        symlink_path.symlink_to(tmpdir, target_is_directory=True)

    # Copy files to tmpdir and decompress them if needed
    if copy_files is not None:
        if hasattr(copy_files, "do_copy"):
            copy_files.do_copy(tmpdir)
        elif isinstance(copy_files, dict):
            for k, v in copy_files.items():
                copy_decompress_files(k, v, tmpdir)

    return tmpdir, job_results_dir


def calc_cleanup(
    atoms: Atoms | None, tmpdir: Path | str, job_results_dir: Path | str
) -> None:
    """
    Perform cleanup operations for a calculation, including gzipping files, copying
    files back to the original directory, and removing the tmpdir.

    Parameters
    ----------
    atoms
        The Atoms object after the calculation. Must have a calculator
        attached. If None, no modifications to the calculator's directory will be made.
    tmpdir
        The path to the tmpdir, where the calculation will be run. It will be
        deleted after the calculation is complete.
    job_results_dir
        The path to the job_results_dir, where the files will ultimately be
        stored.

    Returns
    -------
    None
    """
    job_results_dir, tmpdir = Path(job_results_dir), Path(tmpdir)
    settings = get_settings()

    # Safety check
    tmpdir_base = (settings.SCRATCH_DIR or settings.RESULTS_DIR).resolve()
    if not (
        tmpdir.is_relative_to(tmpdir_base)
        and tmpdir.relative_to(tmpdir_base).parts[0].startswith("tmp-")
    ):
        msg = f"{tmpdir} does not appear to be a tmpdir... exiting for safety!"
        raise ValueError(msg)

    # Update the calculator's directory
    if atoms is not None:
        atoms.calc.directory = job_results_dir

    # Gzip files in tmpdir
    if settings.GZIP_FILES:
        gzip_dir(tmpdir)

    # Move files from tmpdir to job_results_dir.
    LOGGER.info(f"Moving {tmpdir} contents to {job_results_dir}")
    job_results_dir.mkdir(parents=True, exist_ok=True)
    for file_name in os.listdir(tmpdir):
        move(tmpdir / file_name, job_results_dir / file_name)
    rmtree(tmpdir)
    LOGGER.info(f"Calculation results stored at {job_results_dir}")

    # Remove symlink to tmpdir
    if os.name != "nt" and settings.SCRATCH_DIR:
        symlink_path = settings.RESULTS_DIR / f"symlink-{tmpdir.name}"
        symlink_path.unlink(missing_ok=True)


def terminate(tmpdir: Path | str, exception: Exception) -> None:
    """
    Terminate a calculation and move files to a failed directory.

    Parameters
    ----------
    tmpdir
        The path to the tmpdir, where the calculation was run.
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
    settings = get_settings()
    tmpdir = Path(tmpdir)
    job_failed_dir = tmpdir.with_name("failed-" + tmpdir.name.replace("tmp-", ""))
    tmpdir.rename(job_failed_dir)

    msg = f"Calculation failed! Files stored at {job_failed_dir}"
    LOGGER.info(msg)

    if os.name != "nt" and settings.SCRATCH_DIR:
        old_symlink_path = settings.RESULTS_DIR / f"symlink-{tmpdir.name}"
        symlink_path = settings.RESULTS_DIR / f"symlink-{job_failed_dir.name}"
        old_symlink_path.unlink(missing_ok=True)
        symlink_path.symlink_to(job_failed_dir, target_is_directory=True)

    raise JobFailure(job_failed_dir, message=msg, parent_error=exception) from exception
