"""Utility functions for running ASE calculators"""
from __future__ import annotations

import os
from shutil import rmtree
from tempfile import mkdtemp

import numpy as np
from ase import Atoms
from ase.io import Trajectory, read
from ase.optimize import FIRE
from ase.optimize.optimize import Optimizer
from ase.vibrations import Vibrations
from monty.os.path import zpath
from monty.shutil import copy_r, gzip_dir

from quacc import SETTINGS
from quacc.util.atoms import copy_atoms
from quacc.util.files import copy_decompress, make_unique_dir


def run_calc(
    atoms: Atoms,
    geom_file: str | None = None,
    create_unique_workdir: bool | None = None,
    scratch_dir: str | None = None,
    gzip: bool | None = None,
    copy_files: list[str] | None = None,
) -> Atoms:
    """
    Run a calculation in a scratch directory and copy the results back to the
    original directory. This can be useful if file I/O is slow in the working
    directory, so long as file transfer speeds are reasonable.

    This is a wrapper around atoms.get_potential_energy(). Note: This
    function does not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    geom_file
        The filename of the log file that contains the output geometry, used
        to update the atoms object's positions and cell after a job. It is better
        to specify this rather than relying on ASE's atoms.get_potential_energy()
        function to update the positions, as this varies between codes.
    create_unique_workdir
        Whether to automatically create a unique working directory for each calculation.
        If None, defaults to SETTINGS.CREATE_UNIQUE_WORKDIR.
    scratch_dir
        Base path where a tmpdir should be made for running the calculation.
        If None, defaults to SETTINGS.SCRATCH_DIR.
    gzip
        Whether to gzip the output files.
        If None, defaults to SETTINGS.GZIP_FILES.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    Atoms
        The updated Atoms object.
    """

    # Set defaults
    create_unique_workdir = (
        SETTINGS.CREATE_UNIQUE_WORKDIR
        if create_unique_workdir is None
        else create_unique_workdir
    )
    scratch_dir = SETTINGS.SCRATCH_DIR if scratch_dir is None else scratch_dir
    gzip = SETTINGS.GZIP_FILES if gzip is None else gzip

    # Perform staging operations
    start_dir = os.getcwd()
    atoms, tmpdir, results_dir = _calc_setup(
        atoms,
        create_unique_workdir=create_unique_workdir,
        copy_files=copy_files,
        scratch_dir=scratch_dir,
    )

    # Run calculation via get_potential_energy()
    atoms.get_potential_energy()

    # Perform cleanup operations
    _calc_cleanup(start_dir, tmpdir, results_dir, gzip=gzip)

    # Most ASE calculators do not update the atoms object in-place with
    # a call to .get_potential_energy(), which is important if an internal
    # optimizer is used. This section is done to ensure that the atoms object
    # is updated with the correct positions and cell if a `geom_file` is provided.
    if geom_file:
        # Note: We have to be careful to make sure we don't lose the
        # converged magnetic moments, if present. That's why we simply
        # update the positions and cell in-place.
        atoms_new = read(os.path.join(results_dir, zpath(geom_file)))
        if isinstance(atoms_new, list):
            atoms_new = atoms_new[-1]

        # Make sure the atom indices didn't get updated somehow (sanity check). If this
        # happens, there is a serious problem.
        if (
            np.array_equal(atoms_new.get_atomic_numbers(), atoms.get_atomic_numbers())
            is False
        ):
            raise ValueError("Atomic numbers do not match between atoms and geom_file.")

        atoms.positions = atoms_new.positions
        atoms.cell = atoms_new.cell

    return atoms


def run_ase_opt(
    atoms: Atoms,
    fmax: float = 0.01,
    max_steps: int = 500,
    optimizer: Optimizer = FIRE,
    optimizer_kwargs: dict | None = None,
    create_unique_workdir: bool | None = None,
    scratch_dir: str | None = None,
    gzip: bool | None = None,
    copy_files: list[str] | None = None,
) -> Optimizer:
    """
    Run an ASE-based optimization in a scratch directory and copy the results
    back to the original directory. This can be useful if file I/O is slow in
    the working directory, so long as file transfer speeds are reasonable.

    This is a wrapper around the optimizers in ASE. Note: This function does
    not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        Optimizer class to use.
    optimizer_kwargs
        Dictionary of kwargs for the optimizer.
    create_unique_workdir
        Whether to automatically create a unique working directory for each calculation.
        If None, defaults to SETTINGS.CREATE_UNIQUE_WORKDIR.
    scratch_dir
        Base path where a tmpdir should be made for running the calculation.
        If None, defaults to SETTINGS.SCRATCH_DIR.
    gzip
        Whether to gzip the output files.
        If None, defaults to SETTINGS.GZIP_FILES.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    Optimizer
        The ASE Optimizer object.
    """

    # Set defaults
    optimizer_kwargs = optimizer_kwargs or {}
    create_unique_workdir = (
        SETTINGS.CREATE_UNIQUE_WORKDIR
        if create_unique_workdir is None
        else create_unique_workdir
    )
    scratch_dir = SETTINGS.SCRATCH_DIR if scratch_dir is None else scratch_dir
    gzip = SETTINGS.GZIP_FILES if gzip is None else gzip
    start_dir = os.getcwd()

    # Perform staging operations
    atoms, tmpdir, results_dir = _calc_setup(
        atoms,
        create_unique_workdir=create_unique_workdir,
        copy_files=copy_files,
        scratch_dir=scratch_dir,
    )

    # Set Sella kwargs
    if (
        optimizer.__name__ == "Sella"
        and not atoms.pbc.any()
        and "internal" not in optimizer_kwargs
    ):
        optimizer_kwargs["internal"] = True

    # Set up trajectory
    if "trajectory" in optimizer_kwargs:
        raise ValueError("Quacc does not support setting the `trajectory` kwarg.")

    traj_filename = "opt.traj"
    optimizer_kwargs["trajectory"] = Trajectory(traj_filename, "w", atoms=atoms)

    # Define optimizer class
    dyn = optimizer(atoms, **optimizer_kwargs)

    # Run calculation
    dyn.run(fmax=fmax, steps=max_steps)

    # Store the trajectory atoms
    dyn.traj_atoms = read(traj_filename, index=":")

    # Perform cleanup operations
    _calc_cleanup(start_dir, tmpdir, results_dir, gzip=gzip)

    return dyn


def run_ase_vib(
    atoms: Atoms,
    vib_kwargs: dict | None = None,
    create_unique_workdir: bool | None = None,
    scratch_dir: str | None = None,
    gzip: bool | None = None,
    copy_files: list[str] | None = None,
) -> Vibrations:
    """
    Run an ASE-based vibration analysis in a scratch directory and copy the results
    back to the original directory. This can be useful if file I/O is slow in
    the working directory, so long as file transfer speeds are reasonable.

    This is a wrapper around the vibrations module in ASE. Note: This function does
    not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    vib_kwargs
        Dictionary of kwargs for the vibration analysis.
    create_unique_workdir
        Whether to automatically create a unique working directory for each calculation.
        If None, defaults to SETTINGS.CREATE_UNIQUE_WORKDIR.
    scratch_dir
        Base path where a tmpdir should be made for running the calculation.
        If None, defaults to SETTINGS.SCRATCH_DIR.
    gzip
        Whether to gzip the output files.
        If None, defaults to SETTINGS.GZIP_FILES.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    Vibrations
        The updated Vibrations module
    """

    # Set defaults
    vib_kwargs = vib_kwargs or {}
    create_unique_workdir = (
        SETTINGS.CREATE_UNIQUE_WORKDIR
        if create_unique_workdir is None
        else create_unique_workdir
    )
    scratch_dir = SETTINGS.SCRATCH_DIR if scratch_dir is None else scratch_dir
    gzip = SETTINGS.GZIP_FILES if gzip is None else gzip
    start_dir = os.getcwd()

    # Perform staging operations
    atoms, tmpdir, results_dir = _calc_setup(
        atoms,
        create_unique_workdir=create_unique_workdir,
        copy_files=copy_files,
        scratch_dir=scratch_dir,
    )

    # Run calculation
    vib = Vibrations(atoms, name="vib", **vib_kwargs)
    vib.run()
    vib.summary(log="vib_summary.log")

    # Perform cleanup operations
    _calc_cleanup(start_dir, tmpdir, results_dir, gzip=gzip)

    return vib


def _calc_setup(
    atoms: Atoms,
    create_unique_workdir: bool = False,
    copy_files: list[str] | None = None,
    scratch_dir: str | None = None,
) -> tuple[Atoms, str, str]:
    """
    Perform staging operations for a calculation, including copying files
    to the scratch directory, setting the calculator's directory,
    decompressing files, and creating a symlink to the scratch directory.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on with calculator attached.
    create_unique_workdir
        Whether to automatically create a unique working directory for each calculation.
    copy_files
        Filenames to copy from source to scratch directory.
    scratch_dir
        Base path where a tmpdir should be made for running the calculation.

    Returns
    -------
    Atoms
        Copy of the Atoms object with the calculator's directory set.
    str
        The path to the tmpdir, where the calculation will be run. It will be
        deleted after the calculation is complete.
    str
        The path to the results_dir, where the files will ultimately be stored.
        A symlink to the tmpdir will be made here during the calculation for
        convenience.
    """

    if atoms.calc is None:
        raise ValueError("Atoms object must have attached calculator.")

    # Don't modify the original atoms object
    atoms = copy_atoms(atoms)

    # Set where to store the results
    results_dir = make_unique_dir() if create_unique_workdir else os.getcwd()

    # Set the base scratch directory where the tmpdir will be made
    scratch_dir = scratch_dir or os.getcwd()
    if not os.path.exists(scratch_dir):
        os.makedirs(scratch_dir)

    # Create a tmpdir for the calculation within the scratch_dir
    tmpdir = os.path.abspath(mkdtemp(prefix="quacc-tmp-", dir=scratch_dir))

    # Create a symlink (if not on Windows) to the tmpdir in the results_dir
    symlink = os.path.join(results_dir, f"{os.path.basename(tmpdir)}-symlink")
    if os.name != "nt":
        if os.path.islink(symlink):
            os.unlink(symlink)
        os.symlink(tmpdir, symlink)

    # Copy files to tmpdir and decompress them if needed
    if copy_files:
        copy_decompress(copy_files, tmpdir)

    os.chdir(tmpdir)

    return atoms, tmpdir, results_dir


def _calc_cleanup(
    start_dir: str, tmpdir: str, results_dir: str, gzip: bool = True
) -> None:
    """
    Perform cleanup operations for a calculation, including gzipping files,
    copying files back to the original directory, and removing the tmpdir.

    Parameters
    ----------
    start_dir
        The path to the directory where the calculation was started.
    tmpdir
        The path to the tmpdir, where the calculation will be run. It will be
        deleted after the calculation is complete.
    results_dir
        The path to the results_dir, where the files will ultimately be stored.
        A symlink to the tmpdir will be made here during the calculation for
        convenience.
    gzip
        Whether to gzip the output files.

    Returns
    -------
    None
    """

    # Change back to the original directory
    os.chdir(start_dir)

    # Gzip files in tmpdir
    if gzip:
        gzip_dir(tmpdir)

    # Copy files back to results_dir
    copy_r(tmpdir, results_dir)

    # Remove symlink to tmpdir
    symlink = os.path.join(results_dir, f"{os.path.basename(tmpdir)}-symlink")
    if os.path.islink(symlink):
        os.remove(symlink)

    # Remove the tmpdir
    rmtree(tmpdir)
