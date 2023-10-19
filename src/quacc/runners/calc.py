"""Utility functions for running ASE calculators"""
from __future__ import annotations

import os
from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp
from typing import TYPE_CHECKING

import numpy as np
from ase.filters import FrechetCellFilter
from ase.io import Trajectory, read
from ase.optimize import FIRE
from ase.vibrations import Vibrations
from monty.dev import requires
from monty.os.path import zpath
from monty.shutil import copy_r, gzip_dir

from quacc import SETTINGS
from quacc.atoms.core import copy_atoms
from quacc.utils.files import copy_decompress, make_unique_dir

try:
    from sella import Internals, Sella

except ImportError:
    Sella = None

if TYPE_CHECKING:
    from ase import Atoms
    from ase.optimize.optimize import Optimizer


def run_calc(
    atoms: Atoms, geom_file: str | None = None, copy_files: list[str] | None = None
) -> Atoms:
    """
    Run a calculation in a scratch directory and copy the results back to the
    original directory. This can be useful if file I/O is slow in the working
    directory, so long as file transfer speeds are reasonable.

    This is a wrapper around atoms.get_potential_energy(). Note: This function
    does not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    geom_file
        The filename of the log file that contains the output geometry, used to
        update the atoms object's positions and cell after a job. It is better
        to specify this rather than relying on ASE's
        atoms.get_potential_energy() function to update the positions, as this
        varies between codes.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    Atoms
        The updated Atoms object.
    """

    # Perform staging operations
    atoms, tmpdir, job_results_dir = _calc_setup(atoms, copy_files=copy_files)

    # Run calculation via get_potential_energy()
    atoms.get_potential_energy()

    # Most ASE calculators do not update the atoms object in-place with a call
    # to .get_potential_energy(), which is important if an internal optimizer is
    # used. This section is done to ensure that the atoms object is updated with
    # the correct positions and cell if a `geom_file` is provided.
    if geom_file:
        # Note: We have to be careful to make sure we don't lose the converged
        # magnetic moments, if present. That's why we simply update the
        # positions and cell in-place.
        atoms_new = read(zpath(tmpdir / geom_file))
        if isinstance(atoms_new, list):
            atoms_new = atoms_new[-1]

        # Make sure the atom indices didn't get updated somehow (sanity check).
        # If this happens, there is a serious problem.
        if (
            np.array_equal(atoms_new.get_atomic_numbers(), atoms.get_atomic_numbers())
            is False
        ):
            raise ValueError("Atomic numbers do not match between atoms and geom_file.")

        atoms.positions = atoms_new.positions
        atoms.cell = atoms_new.cell

    # Perform cleanup operations
    _calc_cleanup(tmpdir, job_results_dir)

    return atoms


def run_ase_opt(
    atoms: Atoms,
    relax_cell: bool = False,
    fmax: float = 0.01,
    max_steps: int = 500,
    optimizer: Optimizer = FIRE,
    optimizer_kwargs: dict | None = None,
    run_kwargs: dict | None = None,
    copy_files: list[str] | None = None,
) -> Optimizer:
    """
    Run an ASE-based optimization in a scratch directory and copy the results
    back to the original directory. This can be useful if file I/O is slow in
    the working directory, so long as file transfer speeds are reasonable.

    This is a wrapper around the optimizers in ASE. Note: This function does not
    modify the atoms object in-place.

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
        Dictionary of kwargs for the optimizer. Takes all valid kwargs for ASE
        Optimizer classes. Refer to `_set_sella_kwargs` for Sella-related
        kwargs and how they are set.
    run_kwargs
        Dictionary of kwargs for the run() method of the optimizer.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    Optimizer
        The ASE Optimizer object.
    """

    # Set defaults
    optimizer_kwargs = optimizer_kwargs or {}
    run_kwargs = run_kwargs or {}

    # Check if trajectory kwarg is specified
    if "trajectory" in optimizer_kwargs:
        msg = "Quacc does not support setting the `trajectory` kwarg."
        raise ValueError(msg)

    # Perform staging operations
    atoms, tmpdir, job_results_dir = _calc_setup(atoms, copy_files=copy_files)

    # Set Sella kwargs
    if optimizer.__name__ == "Sella":
        _set_sella_kwargs(atoms, optimizer_kwargs)
    optimizer_kwargs.pop("use_TRICs", None)

    # Define the Trajectory object
    traj_filename = tmpdir / "opt.traj"
    traj = Trajectory(traj_filename, "w", atoms=atoms)
    optimizer_kwargs["trajectory"] = traj

    # Set volume relaxation constraints, if relevant
    if relax_cell and atoms.pbc.any():
        atoms = FrechetCellFilter(atoms)

    # Run calculation
    with traj, optimizer(atoms, **optimizer_kwargs) as dyn:
        dyn.run(fmax=fmax, steps=max_steps, **run_kwargs)

    # Store the trajectory atoms
    dyn.traj_atoms = read(traj_filename, index=":")

    # Perform cleanup operations
    _calc_cleanup(tmpdir, job_results_dir)

    return dyn


def run_ase_vib(
    atoms: Atoms, vib_kwargs: dict | None = None, copy_files: list[str] | None = None
) -> Vibrations:
    """
    Run an ASE-based vibration analysis in a scratch directory and copy the
    results back to the original directory. This can be useful if file I/O is
    slow in the working directory, so long as file transfer speeds are
    reasonable.

    This is a wrapper around the vibrations module in ASE. Note: This function
    does not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    vib_kwargs
        Dictionary of kwargs for the vibration analysis.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    Vibrations
        The updated Vibrations module
    """

    # Set defaults
    vib_kwargs = vib_kwargs or {}

    # Perform staging operations
    atoms, tmpdir, job_results_dir = _calc_setup(atoms, copy_files=copy_files)

    # Run calculation
    vib = Vibrations(atoms, name=str(tmpdir / "vib"), **vib_kwargs)
    vib.run()
    vib.summary(log=str(tmpdir / "vib_summary.log"))

    # Perform cleanup operations
    _calc_cleanup(tmpdir, job_results_dir)

    return vib


def _calc_setup(
    atoms: Atoms, copy_files: list[str | Path] | None = None
) -> tuple[Atoms, Path, Path]:
    """
    Perform staging operations for a calculation, including copying files to the
    scratch directory, setting the calculator's directory, decompressing files,
    and creating a symlink to the scratch directory.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on with calculator attached.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    Atoms
        Copy of the Atoms object with the calculator's directory set.
    Path
        The path to the tmpdir, where the calculation will be run. It will be
        deleted after the calculation is complete.
    Path
        The path to the results_dir, where the files will ultimately be stored.
        A symlink to the tmpdir will be made here during the calculation for
        convenience.
    """

    if atoms.calc is None:
        msg = "Atoms object must have attached calculator."
        raise ValueError(msg)

    # Don't modify the original atoms object
    atoms = copy_atoms(atoms)

    # Set where to store the results
    job_results_dir = (
        make_unique_dir(base_path=SETTINGS.RESULTS_DIR)
        if SETTINGS.CREATE_UNIQUE_WORKDIR
        else SETTINGS.RESULTS_DIR
    )

    # Create a tmpdir for the calculation within the scratch_dir
    tmpdir = Path(mkdtemp(prefix="quacc-tmp-", dir=SETTINGS.SCRATCH_DIR)).resolve()

    # Create a symlink to the tmpdir in the results_dir
    if os.name != "nt" and SETTINGS.SCRATCH_DIR != SETTINGS.RESULTS_DIR:
        symlink = job_results_dir / f"{tmpdir.name}-symlink"
        symlink.unlink(missing_ok=True)
        symlink.symlink_to(tmpdir, target_is_directory=True)

    # Copy files to tmpdir and decompress them if needed
    if copy_files:
        copy_decompress(copy_files, tmpdir)

    os.chdir(tmpdir)

    return atoms, tmpdir, job_results_dir


def _calc_cleanup(tmpdir: str | Path, job_results_dir: str | Path) -> None:
    """
    Perform cleanup operations for a calculation, including gzipping files,
    copying files back to the original directory, and removing the tmpdir.

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


@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def _set_sella_kwargs(atoms: Atoms, optimizer_kwargs: dict) -> None:
    """
    Modifies the `optimizer_kwargs` in-place to address various Sella-related
    parameters. This function does the following for the specified key/value
    pairs in `optimizer_kwargs`:

    1. Sets `order = 0` if not specified (i.e. minimization rather than TS
    by default).

    2. If `internal` is not defined and not `atoms.pbc.any()`, set it to `True`.

    3. If `use_TRICs = True` and not `atoms.pbc.any()`, then `internal` is
    built for the user via `find_all_bonds()`, `find_all_angles()`, and
    `find_all_dihedral()`, unless the user has directly specified `internal`.

    Parameters
    ----------
    atoms
        The Atoms object.
    optimizer_kwargs
        The kwargs for the Sella optimizer.

    Returns
    -------
    None
    """

    if "order" not in optimizer_kwargs:
        optimizer_kwargs["order"] = 0

    if not atoms.pbc.any():
        if "internal" not in optimizer_kwargs:
            optimizer_kwargs["internal"] = True

        if optimizer_kwargs.get("use_TRICs") and not isinstance(
            optimizer_kwargs.get("internal"), Internals
        ):
            internals = Internals(atoms, allow_fragments=True)
            internals.find_all_bonds()
            internals.find_all_angles()
            internals.find_all_dihedrals()
            optimizer_kwargs["internal"] = internals
