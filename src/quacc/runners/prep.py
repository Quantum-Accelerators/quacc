"""Prepration for runners"""
from __future__ import annotations

import os
from datetime import datetime, timezone
from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp

from monty.shutil import copy_r, gzip_dir

from quacc import SETTINGS
from quacc.utils.files import copy_decompress, make_unique_dir


def calc_setup(copy_files: list[str | Path] | None = None) -> tuple[Path, Path]:
    """
    Perform staging operations for a calculation, including copying files to the scratch
    directory, setting the calculator's directory, decompressing files, and creating a
    symlink to the scratch directory.

    Parameters
    ----------
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    Path
        The path to the tmpdir, where the calculation will be run. It will be
        deleted after the calculation is complete.
    Path
        The path to the results_dir, where the files will ultimately be stored.
        A symlink to the tmpdir will be made here during the calculation for
        convenience.
    """

    # Set where to store the results
    job_results_dir = (
        make_unique_dir(base_path=SETTINGS.RESULTS_DIR)
        if SETTINGS.CREATE_UNIQUE_WORKDIR
        else SETTINGS.RESULTS_DIR
    )

    # Create a tmpdir for the calculation within the scratch_dir
    time_now = datetime.now(timezone.utc).strftime("%Y-%m-%d-%H-%M-%S-%f")
    tmpdir = Path(
        mkdtemp(prefix=f"quacc-tmp-{time_now}-", dir=SETTINGS.SCRATCH_DIR)
    ).resolve()

    # Create a symlink to the tmpdir in the results_dir
    if os.name != "nt" and SETTINGS.SCRATCH_DIR != SETTINGS.RESULTS_DIR:
        symlink = job_results_dir / f"{tmpdir.name}-symlink"
        symlink.unlink(missing_ok=True)
        symlink.symlink_to(tmpdir, target_is_directory=True)

    # Copy files to tmpdir and decompress them if needed
    if copy_files:
        copy_decompress(copy_files, tmpdir)

    os.chdir(tmpdir)

    return tmpdir, job_results_dir


def calc_cleanup(tmpdir: str | Path, job_results_dir: str | Path) -> None:
    """
    Perform cleanup operations for a calculation, including gzipping files, copying
    files back to the original directory, and removing the tmpdir.

    Parameters
    ----------
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

    # Change to the results directory
    os.chdir(job_results_dir)

    # Gzip files in tmpdir
    if SETTINGS.GZIP_FILES:
        gzip_dir(tmpdir)

    # Copy files back to job_results_dir
    copy_r(tmpdir, job_results_dir)

    # Remove symlink to tmpdir
    symlink_path = job_results_dir / f"{tmpdir.name}-symlink"
    symlink_path.unlink(missing_ok=True)

    # Remove the tmpdir
    rmtree(tmpdir, ignore_errors=True)
