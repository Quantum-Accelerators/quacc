"""Utility functions for running ASE calculators with ASE-based methods."""

from __future__ import annotations

import sys
from shutil import copy, copytree
from typing import TYPE_CHECKING

import numpy as np
from ase.filters import FrechetCellFilter
from ase.io import Trajectory, read
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from monty.dev import requires
from monty.os.path import zpath

from quacc import SETTINGS
from quacc.atoms.core import copy_atoms, get_final_atoms_from_dyn
from quacc.runners.prep import calc_cleanup, calc_setup
from quacc.utils.dicts import recursive_dict_merge

try:
    from sella import Sella

except ImportError:
    Sella = None

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any, TypedDict

    from ase.atoms import Atoms
    from ase.optimize.optimize import Optimizer

    from quacc.utils.files import Filenames, SourceDirectory

    class OptimizerKwargs(TypedDict, total=False):
        """
        Type hint for `optimizer_kwargs` in [quacc.runners.ase.run_opt][].
        """

        restart: Path | str | None  # default = None
        append_trajectory: bool  # default = False

    class VibKwargs(TypedDict, total=False):
        """
        Type hint for `vib_kwargs` in [quacc.runners.ase.run_vib][].
        """

        indices: list[int] | None  # default = None
        delta: float  # default = 0.01
        nfree: int  # default = 2


def run_calc(
    atoms: Atoms,
    geom_file: str | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    get_forces: bool = False,
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
        Files to copy (and decompress) from source to the runtime directory.
    get_forces
        Whether to use `atoms.get_forces()` instead of `atoms.get_potential_energy()`.

    Returns
    -------
    Atoms
        The updated Atoms object.
    """

    # Copy atoms so we don't modify it in-place
    atoms = copy_atoms(atoms)

    # Perform staging operations
    tmpdir, job_results_dir = calc_setup(atoms, copy_files=copy_files)

    # Run calculation
    if get_forces:
        atoms.get_forces()
    else:
        atoms.get_potential_energy()

    # Most ASE calculators do not update the atoms object in-place with a call
    # to .get_potential_energy(), which is important if an internal optimizer is
    # used. This section is done to ensure that the atoms object is updated to
    # the final geometry if `geom_file` is provided.
    # Note: We have to be careful to make sure we don't lose the calculator
    # object, as this contains important information such as the parameters
    # and output properties (e.g. final magnetic moments).
    if geom_file:
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
    calc_cleanup(atoms, tmpdir, job_results_dir)

    return atoms


def run_opt(
    atoms: Atoms,
    relax_cell: bool = False,
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: Optimizer = BFGS,
    optimizer_kwargs: OptimizerKwargs | None = None,
    store_intermediate_files: bool = False,
    run_kwargs: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
    relax_cell
        Whether to relax the unit cell shape and volume.
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
    store_intermediate_files
        Whether to store the files generated at each intermediate step in the
        optimization. If enabled, they will be stored in a directory named
        `stepN` where `N` is the step number, starting at 0.
    run_kwargs
        Dictionary of kwargs for the run() method of the optimizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    Optimizer
        The ASE Optimizer object.
    """

    # Copy atoms so we don't modify it in-place
    atoms = copy_atoms(atoms)

    # Perform staging operations
    tmpdir, job_results_dir = calc_setup(atoms, copy_files=copy_files)

    # Set defaults
    optimizer_kwargs = recursive_dict_merge(
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

    # Handle optimizer kwargs
    if optimizer.__name__.startswith("SciPy"):
        optimizer_kwargs.pop("restart")
    elif optimizer.__name__ == "Sella":
        _set_sella_kwargs(atoms, optimizer_kwargs)
    elif optimizer.__name__ == "IRC":
        optimizer_kwargs.pop("restart", None)

    # Define the Trajectory object
    traj_file = tmpdir / "opt.traj"
    traj = Trajectory(traj_file, "w", atoms=atoms)
    optimizer_kwargs["trajectory"] = traj

    # Set volume relaxation constraints, if relevant
    if relax_cell and atoms.pbc.any():
        atoms = FrechetCellFilter(atoms)

    # Run optimization
    with traj, optimizer(atoms, **optimizer_kwargs) as dyn:
        if store_intermediate_files:
            opt = dyn.irun(fmax=fmax, steps=max_steps, **run_kwargs)
            for i, _ in enumerate(opt):
                _copy_intermediate_files(
                    tmpdir,
                    i,
                    files_to_ignore=[
                        traj_file,
                        optimizer_kwargs["restart"],
                        optimizer_kwargs["logfile"],
                    ],
                )
        else:
            dyn.run(fmax=fmax, steps=max_steps, **run_kwargs)

    # Store the trajectory atoms
    dyn.traj_atoms = read(traj_file, index=":")

    # Perform cleanup operations
    calc_cleanup(get_final_atoms_from_dyn(dyn), tmpdir, job_results_dir)

    return dyn


def run_vib(
    atoms: Atoms,
    vib_kwargs: VibKwargs | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        Dictionary of kwargs for the [ase.vibrations.Vibrations][] class.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

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
    tmpdir, job_results_dir = calc_setup(atoms, copy_files=copy_files)

    # Run calculation
    vib = Vibrations(atoms, name=str(tmpdir / "vib"), **vib_kwargs)
    vib.run()

    # Summarize run
    vib.summary(log=sys.stdout if SETTINGS.DEBUG else str(tmpdir / "vib_summary.log"))

    # Perform cleanup operations
    calc_cleanup(vib.atoms, tmpdir, job_results_dir)

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

    if not atoms.pbc.any() and "internal" not in optimizer_kwargs:
        optimizer_kwargs["internal"] = True


def _copy_intermediate_files(
    tmpdir: Path, step_number: int, files_to_ignore: list[Path] | None = None
) -> None:
    """
    Copy all files in the working directory to a subdirectory named `stepN` where `N`
    is the step number. This is useful for storing intermediate files generated during
    an ASE relaaxation.

    Parameters
    ----------
    tmpdir
        The working directory.
    step_number
        The step number.
    files_to_ignore
        A list of files to ignore when copying files to the subdirectory.

    Returns
    -------
    None
    """
    files_to_ignore = files_to_ignore or []
    store_path = tmpdir / f"step{step_number}"
    store_path.mkdir()
    for item in tmpdir.iterdir():
        if not item.name.startswith("step") and item not in files_to_ignore:
            if item.is_file():
                copy(item, store_path)
            elif item.is_dir():
                copytree(item, store_path / item.name)
