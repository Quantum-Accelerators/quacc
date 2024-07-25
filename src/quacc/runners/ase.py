"""Utility functions for running ASE calculators with ASE-based methods."""

from __future__ import annotations

import logging
import sys
from importlib.util import find_spec
from shutil import copy, copytree
from typing import TYPE_CHECKING, Any

import numpy as np
from ase.atoms import Atoms
from ase.calculators import calculator
from ase.filters import FrechetCellFilter
from ase.io import Trajectory, read
from ase.md.md import MolecularDynamics
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation,
)
from ase.mep import NEB
from ase.mep.neb import NEBOptimizer
from ase.optimize import BFGS
from ase.optimize.sciopt import SciPyOptimizer
from ase.vibrations import Vibrations
from monty.dev import requires
from monty.os.path import zpath

from quacc import get_settings
from quacc.atoms.core import copy_atoms
from quacc.runners._base import BaseRunner
from quacc.runners.prep import calc_cleanup, calc_setup, terminate
from quacc.utils.dicts import recursive_dict_merge

LOGGER = logging.getLogger(__name__)

has_sella = bool(find_spec("sella"))
has_geodesic_interpolate = bool(find_spec("geodesic_interpolate"))

if has_geodesic_interpolate:
    from geodesic_interpolate.geodesic import Geodesic
    from geodesic_interpolate.interpolation import redistribute

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any, Literal

    from ase.calculators.calculator import Calculator
    from ase.optimize.optimize import Dynamics, Optimizer

    from quacc.types import (
        Filenames,
        MaxwellBoltzmanDistributionKwargs,
        SourceDirectory,
        VibKwargs,
    )


class Runner(BaseRunner):
    """
    Run various types of calculations in a scratch directory and copy the results back
    to the original directory. Note: This function does not modify the atoms object in-place.
    """

    def __init__(
        self,
        atoms: Atoms,
        calculator: Calculator,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    ) -> None:
        """
        Initialize the Runner object.

        Parameters
        ----------
        atoms
            The Atoms object to run calculations on.
        calculator
            The instantiated ASE calculator object to attach to the Atoms object.
        copy_files
            Files to copy (and decompress) from source to the runtime directory.

        Returns
        -------
        None
        """
        self.atoms = copy_atoms(atoms)
        self.atoms.calc = calculator
        self.copy_files = copy_files
        self.setup()

    def run_calc(
        self, geom_file: str | None = None, properties: list[str] | None = None
    ) -> Atoms:
        """
        This is a wrapper around `atoms.calc.calculate()`.

        Parameters
        ----------
        geom_file
            The filename of the log file that contains the output geometry, used to
            update the atoms object's positions and cell after a job. It is better
            to specify this rather than relying on ASE to update the positions, as the
            latter behavior varies between codes.
        properties
            List of properties to calculate. Defaults to ["energy"] if `None`.

        Returns
        -------
        Atoms
            The updated Atoms object.
        """
        if properties is None:
            properties = ["energy"]

        # Run calculation
        try:
            self.atoms.calc.calculate(self.atoms, properties, calculator.all_changes)
        except Exception as exception:
            terminate(self.tmpdir, exception)

        # Most ASE calculators do not update the atoms object in-place with a call
        # to .get_potential_energy(), which is important if an internal optimizer is
        # used. This section is done to ensure that the atoms object is updated to
        # the final geometry if `geom_file` is provided.
        # Note: We have to be careful to make sure we don't lose the calculator
        # object, as this contains important information such as the parameters
        # and output properties (e.g. final magnetic moments).
        if geom_file:
            atoms_new = read(zpath(str(self.tmpdir / geom_file)))
            if isinstance(atoms_new, list):
                atoms_new = atoms_new[-1]

            # Make sure the atom indices didn't get updated somehow (sanity check).
            # If this happens, there is a serious problem.
            if (
                np.array_equal(
                    atoms_new.get_atomic_numbers(), self.atoms.get_atomic_numbers()
                )
                is False
            ):
                raise ValueError(
                    "Atomic numbers do not match between atoms and geom_file."
                )

            self.atoms.positions = atoms_new.positions
            self.atoms.cell = atoms_new.cell

        # Perform cleanup operations
        self.cleanup()

        return self.atoms

    def run_opt(
        self,
        relax_cell: bool = False,
        fmax: float | None = 0.01,
        max_steps: int = 1000,
        optimizer: Dynamics = BFGS,
        optimizer_kwargs: dict[str, Any] | None = None,
        store_intermediate_results: bool = False,
        fn_hook: callable | None = None,
        run_kwargs: dict[str, Any] | None = None,
    ) -> Dynamics:
        """
        This is a wrapper around the optimizers in ASE.

        Parameters
        ----------
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
        store_intermediate_results
            Whether to store the files generated at each intermediate step in the
            optimization. If enabled, they will be stored in a directory named
            `stepN` where `N` is the step number, starting at 0.
        fn_hook
            A custom function to call after each step of the optimization.
            The function must take the instantiated dynamics class as
            its only argument.
        run_kwargs
            Dictionary of kwargs for the `run()` method of the optimizer.

        Returns
        -------
        Dynamics
            The ASE Dynamics object following an optimization.
        """
        # Set defaults
        settings = get_settings()
        merged_optimizer_kwargs = recursive_dict_merge(
            {
                "logfile": "-" if settings.DEBUG else self.tmpdir / "opt.log",
                "restart": self.tmpdir / "opt.json",
            },
            optimizer_kwargs,
        )
        if merged_optimizer_kwargs is None:
            merged_optimizer_kwargs = {}
        run_kwargs = run_kwargs or {}
        traj_filename = "opt.traj"

        # Check if trajectory kwarg is specified
        if "trajectory" in merged_optimizer_kwargs:
            msg = "Quacc does not support setting the `trajectory` kwarg."
            raise ValueError(msg)

        # Handle optimizer kwargs
        if (
            issubclass(optimizer, (SciPyOptimizer, MolecularDynamics))
            or optimizer.__name__ == "IRC"
        ):
            # https://gitlab.com/ase/ase/-/issues/1476
            # https://gitlab.com/ase/ase/-/merge_requests/3310
            merged_optimizer_kwargs.pop("restart", None)
        if optimizer.__name__ == "Sella":
            self._set_sella_kwargs(merged_optimizer_kwargs)

        # Define the Trajectory object
        traj_file = self.tmpdir / traj_filename
        traj = Trajectory(traj_file, "w", atoms=self.atoms)
        merged_optimizer_kwargs["trajectory"] = traj

        # Set volume relaxation constraints, if relevant
        if relax_cell and self.atoms.pbc.any():
            self.atoms = FrechetCellFilter(self.atoms)

        # Run optimization
        full_run_kwargs = {"fmax": fmax, "steps": max_steps, **run_kwargs}
        if issubclass(optimizer, MolecularDynamics):
            full_run_kwargs.pop("fmax")
        try:
            with traj, optimizer(self.atoms, **merged_optimizer_kwargs) as dyn:
                if issubclass(optimizer, (SciPyOptimizer, MolecularDynamics)):
                    # https://gitlab.coms/ase/ase/-/issues/1475
                    # https://gitlab.com/ase/ase/-/issues/1497
                    dyn.run(**full_run_kwargs)
                else:
                    for i, _ in enumerate(dyn.irun(**full_run_kwargs)):
                        if store_intermediate_results:
                            self._copy_intermediate_files(
                                i,
                                files_to_ignore=[
                                    traj_file,
                                    merged_optimizer_kwargs.get("restart"),
                                    merged_optimizer_kwargs.get("logfile"),
                                ],
                            )
                        if fn_hook:
                            fn_hook(dyn)
        except Exception as exception:
            terminate(self.tmpdir, exception)

        # Perform cleanup operations
        self.cleanup()
        traj.filename = zpath(str(self.job_results_dir / traj_filename))
        dyn.trajectory = traj

        return dyn

    def run_vib(self, vib_kwargs: VibKwargs | None = None) -> Vibrations:
        """
        Run an ASE-based vibration analysis in a scratch directory and copy the results back
        to the original directory. This can be useful if file I/O is slow in the working
        directory, so long as file transfer speeds are reasonable.

        This is a wrapper around the vibrations module in ASE.

        Parameters
        ----------
        vib_kwargs
            Dictionary of kwargs for the [ase.vibrations.Vibrations][] class.

        Returns
        -------
        Vibrations
            The updated Vibrations module
        """
        # Set defaults
        vib_kwargs = vib_kwargs or {}
        settings = get_settings()

        # Run calculation
        vib = Vibrations(self.atoms, name=str(self.tmpdir / "vib"), **vib_kwargs)
        try:
            vib.run()
        except Exception as exception:
            terminate(self.tmpdir, exception)

        # Summarize run
        vib.summary(
            log=sys.stdout if settings.DEBUG else str(self.tmpdir / "vib_summary.log")
        )

        # Perform cleanup operations
        self.cleanup()

        return vib

    def run_md(
        self,
        dynamics: MolecularDynamics,
        dynamics_kwargs: dict[str, Any] | None = None,
        steps: int = 1000,
        maxwell_boltzmann_kwargs: MaxwellBoltzmanDistributionKwargs | None = None,
        set_com_stationary: bool = False,
        set_zero_rotation: bool = False,
    ) -> MolecularDynamics:
        """
        Run an ASE-based MD in a scratch directory and copy the results back to
        the original directory.

        Parameters
        ----------
        dynamics
            MolecularDynamics class to use, from `ase.md.md.MolecularDynamics`.
        dynamics_kwargs
            Dictionary of kwargs for the dynamics. Takes all valid kwargs for ASE
            MolecularDynamics classes.
        steps
            Maximum number of steps to run
        maxwell_boltzmann_kwargs
            If specified, a `MaxwellBoltzmannDistribution` will be applied to the atoms
            based on `ase.md.velocitydistribution.MaxwellBoltzmannDistribution` with the
            specified keyword arguments.
        set_com_stationary
            Whether to set the center-of-mass momentum to zero. This would be applied after
            any `MaxwellBoltzmannDistribution` is set.
        set_zero_rotation
            Whether to set the total angular momentum to zero. This would be applied after
            any `MaxwellBoltzmannDistribution` is set.

        Returns
        -------
        MolecularDymamics
            The ASE MolecularDynamics object.
        """

        # Set defaults
        dynamics_kwargs = dynamics_kwargs or {}
        maxwell_boltzmann_kwargs = maxwell_boltzmann_kwargs or {}
        settings = get_settings()
        dynamics_kwargs["logfile"] = "-" if settings.DEBUG else self.tmpdir / "md.log"

        if maxwell_boltzmann_kwargs:
            MaxwellBoltzmannDistribution(self.atoms, **maxwell_boltzmann_kwargs)
        if set_com_stationary:
            Stationary(self.atoms)
        if set_zero_rotation:
            ZeroRotation(self.atoms)

        return self.run_opt(
            fmax=None,
            max_steps=steps,
            optimizer=dynamics,
            optimizer_kwargs=dynamics_kwargs,
        )

    def _copy_intermediate_files(
        self, step_number: int, files_to_ignore: list[Path] | None = None
    ) -> None:
        """
        Copy all files in the working directory to a subdirectory named `stepN` where `N`
        is the step number. This is useful for storing intermediate files generated during
        an ASE relaaxation.

        Parameters
        ----------
        step_number
            The step number.
        files_to_ignore
            A list of files to ignore when copying files to the subdirectory.

        Returns
        -------
        None
        """
        files_to_ignore = files_to_ignore or []
        store_path = self.tmpdir / f"step{step_number}"
        store_path.mkdir()
        for item in self.tmpdir.iterdir():
            if not item.name.startswith("step") and item not in files_to_ignore:
                if item.is_file():
                    copy(item, store_path)
                elif item.is_dir():
                    copytree(item, store_path / item.name)

    @requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
    def _set_sella_kwargs(self, optimizer_kwargs: dict[str, Any]) -> None:
        """
        Modifies the `optimizer_kwargs` in-place to address various Sella-related
        parameters. This function does the following for the specified key/value pairs in
        `optimizer_kwargs`:

        1. Sets `order = 0` if not specified (i.e. minimization rather than TS
        by default).

        2. If `internal` is not defined and not `atoms.pbc.any()`, set it to `True`.

        Parameters
        ----------
        optimizer_kwargs
            The kwargs for the Sella optimizer.

        Returns
        -------
        None
        """
        if "order" not in optimizer_kwargs:
            optimizer_kwargs["order"] = 0

        if not self.atoms.pbc.any() and "internal" not in optimizer_kwargs:
            optimizer_kwargs["internal"] = True


def run_neb(
    images: list[Atoms],
    relax_cell: bool = False,
    fmax: float = 0.01,
    max_steps: int | None = 1000,
    optimizer: NEBOptimizer | Optimizer = NEBOptimizer,
    optimizer_kwargs: dict[str, Any] | None = None,
    neb_kwargs: dict[str, Any] | None = None,
    run_kwargs: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> Dynamics:
    """
    Run NEB optimization.

    Parameters
    ----------
    images
        List of images representing the initial path.
    relax_cell
        Whether to relax the unit cell shape and volume.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        Optimizer class to use. All Optimizers except BFGSLineSearch
    optimizer_kwargs
        Dictionary of kwargs for the optimizer.
    run_kwargs
        Dictionary of kwargs for the run() method of the optimizer.
    neb_kwargs
        Dictionary of kwargs for the NEB.
    copy_files
        Files to copy before running the calculation.

    Returns
    -------
    Dynamics
        The ASE Dynamics object following an optimization.
    """
    if optimizer.__name__ == "BFGSLineSearch":
        raise ValueError("BFGSLineSearch is not allowed as optimizer with NEB.")

    # Copy atoms so we don't modify it in-place
    images = copy_atoms(images)
    settings = get_settings()

    neb = NEB(images, **neb_kwargs)

    dir_lists = []
    # Perform staging operations
    # this calc_setup function is not suited for multiple Atoms objects
    for image in images:
        tmpdir_i, job_results_dir_i = calc_setup(image, copy_files=copy_files)
        dir_lists.append([tmpdir_i, job_results_dir_i])

    # Set defaults
    optimizer_kwargs = recursive_dict_merge(
        {
            "logfile": "-" if settings.DEBUG else dir_lists[0][1] / "opt.log",
            "restart": dir_lists[0][1] / "opt.json",
        },
        optimizer_kwargs,
    )
    run_kwargs = run_kwargs or {}
    traj_filename = "opt.traj"
    # Check if trajectory kwarg is specified
    if "trajectory" in optimizer_kwargs:
        msg = "Quacc does not support setting the `trajectory` kwarg."
        raise ValueError(msg)

    # Define the Trajectory object
    traj_file = dir_lists[0][0] / traj_filename
    traj = Trajectory(traj_file, "w", atoms=neb)

    # Set volume relaxation constraints, if relevant
    if relax_cell:
        for i in range(len(images)):
            if images[i].pbc.any():
                images[i] = FrechetCellFilter(images[i])

    dyn = optimizer(neb, **optimizer_kwargs)
    dyn.attach(traj.write)
    dyn.run(fmax, max_steps)
    traj.close()

    # Perform cleanup operations first images's results directory contains traj file.
    for i, image in enumerate(images):
        calc_cleanup(image, dir_lists[i][0], dir_lists[i][1])

    traj.filename = str(zpath(dir_lists[0][1] / traj_filename))
    dyn.trajectory = traj
    return dyn


@requires(
    has_geodesic_interpolate,
    "geodesic-interpolate must be installed. Refer to the quacc documentation.",
)
def _geodesic_interpolate_wrapper(
    reactant: Atoms,
    product: Atoms,
    n_images: int = 20,
    perform_sweep: bool | Literal["auto"] = "auto",
    convergence_tolerance1: float = 1e-2,
    convergence_tolerance2: float = 2e-3,
    max_iterations: int = 15,
    max_micro_iterations: int = 20,
    morse_scaling: float = 1.7,
    geometry_friction: float = 1e-2,
    distance_cutoff: float = 3.0,
    sweep_cutoff_size: int = 35,
) -> list[Atoms]:
    """
    Interpolates between two geometries and optimizes the path with the geodesic method.

    Parameters
    ----------
    reactant
        The ASE Atoms object representing the initial geometry.
    product
        The ASE Atoms object representing the final geometry.
    n_images
        Number of images for interpolation. Default is 20.
    perform_sweep
        Whether to sweep across the path optimizing one image at a time.
        Default is to perform sweeping updates if there are more than 35 atoms.
    convergence_tolerance1
        the value passed to the tol keyword argument of
         geodesic_interpolate.interpolation.redistribute. Default is 1e-2.
    convergence_tolerance2
        the value passed to the tol keyword argument of geodesic_smoother.smooth
        or geodesic_smoother.sweep. Default is 2e-3.
    max_iterations
        Maximum number of minimization iterations. Default is 15.
    max_micro_iterations
        Maximum number of micro iterations for the sweeping algorithm. Default is 20.
    morse_scaling
        Exponential parameter for the Morse potential. Default is 1.7.
    geometry_friction
        Size of friction term used to prevent very large changes in geometry. Default is 1e-2.
    distance_cutoff
        Cut-off value for the distance between a pair of atoms to be included in the coordinate system. Default is 3.0.
    sweep_cutoff_size
        Cut off system size that above which sweep function will be called instead of smooth
        in Geodesic.

    Returns
    -------
    list[Atoms]
        A list of ASE Atoms objects representing the smoothed path between the reactant and product geometries.
    """
    reactant = copy_atoms(reactant)
    product = copy_atoms(product)

    # Read the initial geometries.
    chemical_symbols = reactant.get_chemical_symbols()

    # First redistribute number of images. Perform interpolation if too few and subsampling if too many images are given
    raw_interpolated_positions = redistribute(
        chemical_symbols,
        [reactant.positions, product.positions],
        n_images,
        tol=convergence_tolerance1,
    )

    # Perform smoothing by minimizing distance in Cartesian coordinates with redundant internal metric
    # to find the appropriate geodesic curve on the hyperspace.
    geodesic_smoother = Geodesic(
        chemical_symbols,
        raw_interpolated_positions,
        morse_scaling,
        threshold=distance_cutoff,
        friction=geometry_friction,
    )
    if perform_sweep == "auto":
        perform_sweep = len(chemical_symbols) > sweep_cutoff_size
    if perform_sweep:
        geodesic_smoother.sweep(
            tol=convergence_tolerance2,
            max_iter=max_iterations,
            micro_iter=max_micro_iterations,
        )
    else:
        geodesic_smoother.smooth(tol=convergence_tolerance2, max_iter=max_iterations)
    return [
        Atoms(symbols=chemical_symbols, positions=geom)
        for geom in geodesic_smoother.path
    ]
