"""Utility functions for running ASE calculators with ASE-based methods."""
from __future__ import annotations

import sys
from typing import TYPE_CHECKING

import numpy as np
from ase.filters import FrechetCellFilter
from ase.io import Trajectory, read
from ase.optimize import FIRE
from ase.vibrations import Vibrations
from monty.dev import requires
from monty.os.path import zpath

from quacc import SETTINGS
from quacc.atoms.core import copy_atoms
from quacc.runners.prep import calc_cleanup, calc_setup
from quacc.utils.dicts import merge_dicts

try:
    from sella import Internals, Sella

except ImportError:
    Sella = None

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any, TypedDict

    from ase.atoms import Atoms
    from ase.optimize.optimize import Optimizer

    class OptimizerKwargs(TypedDict, total=False):
        restart: str | None  # default = None
        append_trajectory: bool  # default = False

    class VibKwargs(TypedDict, total=False):
        indices: list[int] | None  # default = None
        delta: float  # default = 0.01
        nfree: int  # default = 2


def run_calc(
    atoms: Atoms,
    geom_file: str | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
) -> Atoms:
    """
    Run a calculation in a scratch directory and copy the results back to the original
    directory. This can be useful if file I/O is slow in the working directory, so long
    as file transfer speeds are reasonable.

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

    # Copy atoms so we don't modify it in-place
    atoms = copy_atoms(atoms)

    # Perform staging operations
    tmpdir, job_results_dir = calc_setup(copy_files=copy_files)

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
    calc_cleanup(tmpdir, job_results_dir)

    return atoms


def run_opt(
    atoms: Atoms,
    relax_cell: bool = False,
    fmax: float = 0.01,
    max_steps: int = 500,
    optimizer: Optimizer = FIRE,
    optimizer_kwargs: OptimizerKwargs | None = None,
    run_kwargs: dict[str, Any] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
) -> Optimizer:
    """
    Run an ASE-based optimization in a scratch directory and copy the results back to
    the original directory. This can be useful if file I/O is slow in the working
    directory, so long as file transfer speeds are reasonable.

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

    # Copy atoms so we don't modify it in-place
    atoms = copy_atoms(atoms)

    # Perform staging operations
    tmpdir, job_results_dir = calc_setup(copy_files=copy_files)

    # Set defaults
    optimizer_kwargs = merge_dicts(
        {
            "logfile": "-" if SETTINGS.DEBUG else tmpdir / "opt.log",
            "restart": tmpdir / "opt.pckl",
        },
        optimizer_kwargs,
    )
    run_kwargs = run_kwargs or {}

    # Check if trajectory kwarg is specified
    if "trajectory" in optimizer_kwargs:
        msg = "Quacc does not support setting the `trajectory` kwarg."
        raise ValueError(msg)

    # Set Sella kwargs
    if optimizer.__name__ == "Sella":
        _set_sella_kwargs(atoms, optimizer_kwargs)
    elif optimizer.__name__ == "IRC":
        optimizer_kwargs.pop("restart", None)
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
    calc_cleanup(tmpdir, job_results_dir)

    return dyn


def run_vib(
    atoms: Atoms,
    vib_kwargs: VibKwargs | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
) -> Vibrations:
    """
    Run an ASE-based vibration analysis in a scratch directory and copy the results back
    to the original directory. This can be useful if file I/O is slow in the working
    directory, so long as file transfer speeds are reasonable.

    This is a wrapper around the vibrations module in ASE. Note: This function
    does not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    vib_kwargs
        Dictionary of kwargs for the `ase.vibrations.Vibrations` class.
    copy_files
        Filenames to copy from source to scratch directory.

    Returns
    -------
    Vibrations
        The updated Vibrations module
    """

    # Copy atoms so we don't modify it in-place
    atoms = copy_atoms(atoms)

    # Set defaults
    vib_kwargs = vib_kwargs or {}

    # Perform staging operations
    tmpdir, job_results_dir = calc_setup(copy_files=copy_files)

    # Run calculation
    vib = Vibrations(atoms, name=str(tmpdir / "vib"), **vib_kwargs)
    vib.run()

    # Summarize run
    vib.summary(log=sys.stdout if SETTINGS.DEBUG else str(tmpdir / "vib_summary.log"))

    # Perform cleanup operations
    calc_cleanup(tmpdir, job_results_dir)

    return vib


@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def _set_sella_kwargs(atoms: Atoms, optimizer_kwargs: dict[str, Any]) -> None:
    """
    Modifies the `optimizer_kwargs` in-place to address various Sella-related
    parameters. This function does the following for the specified key/value pairs in
    `optimizer_kwargs`:

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
