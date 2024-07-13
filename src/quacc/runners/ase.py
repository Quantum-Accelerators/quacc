"""Utility functions for running ASE calculators with ASE-based methods."""

from __future__ import annotations

import logging
import sys
from importlib.util import find_spec
from shutil import copy, copytree
from typing import TYPE_CHECKING, Callable

import numpy as np
from ase.calculators import calculator
from ase.filters import FrechetCellFilter
from ase.io import Trajectory, read
from ase.md.md import MolecularDynamics
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation,
)
from ase.optimize import BFGS
from ase.optimize.sciopt import SciPyOptimizer
from ase.vibrations import Vibrations
from monty.dev import requires
from monty.os.path import zpath

from quacc import get_settings
from quacc.atoms.core import copy_atoms
from quacc.runners._base import BaseRunner
from quacc.runners.prep import terminate
from quacc.utils.dicts import recursive_dict_merge

LOGGER = logging.getLogger(__name__)

has_sella = bool(find_spec("sella"))


if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any

    from ase.atoms import Atoms
    from ase.calculators.calculator import Calculator
    from ase.optimize.optimize import Dynamics

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
        fn_hook: Callable | None = None,
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
