"""Core functions for common recipes."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

from ase.optimize import BFGS

from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize, VibSummarize

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from ase.atoms import Atoms
    from ase.calculators.calculator import BaseCalculator
    from ase.optimize.optimize import Optimizer

    from quacc.types import Filenames, OptSchema, RunSchema, SourceDirectory, VibSchema


class Recipe:
    """Base class for all discrete compute tasks."""

    def __init__(self, calculator: BaseCalculator) -> None:
        """Initialize the recipe.

        Parameters
        ----------
        calculator
            The instantiated ASE calculator
        """
        self.calculator = calculator

    def run(
        self,
        atoms: Atoms,
        geom_file: str | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        additional_fields: dict[str, Any] | None = None,
    ) -> RunSchema:
        """Run the executable with the parameters.

        Parameters
        ----------
        atoms
            Atoms object
        geom_file
            The filename of the log file that contains the output geometry, used to
            update the atoms object's positions and cell after a job. It is better
            to specify this rather than relying on ASE to update the positions, as the
            latter behavior varies between codes.
        copy_files
            Files to copy to runtime directory
        additional_fields
            Additional fields for results
        **calc_kwargs
            Calculator parameters that override defaults

        Returns
        -------
        RunSchema
            Results dictionary
        """
        additional_fields = additional_fields or {}
        additional_fields = {
            "name": f"{self.calculator.__class__.__name__}"
        } | additional_fields

        final_atoms = Runner(atoms, self.calculator, copy_files=copy_files).run_calc(
            geom_file=geom_file
        )
        return Summarize(additional_fields=additional_fields).run(final_atoms, atoms)

    def static(
        self,
        atoms: Atoms,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        additional_fields: dict[str, Any] | None = None,
    ) -> RunSchema:
        """Run a static calculation.

        Parameters
        ----------
        atoms
            Atoms object
        copy_files
            Files to copy to runtime directory
        additional_fields
            Additional fields for results

        Returns
        -------
        RunSchema
            Results dictionary
        """
        additional_fields = additional_fields or {}
        additional_fields = {
            "name": f"{self.calculator.__class__.__name__} Static"
        } | additional_fields

        return self.run(
            atoms, copy_files=copy_files, additional_fields=additional_fields
        )

    def relax(
        self,
        atoms: Atoms,
        relax_cell: bool = False,
        fmax: float | None = 0.01,
        max_steps: int = 1000,
        optimizer: type[Optimizer] = BFGS,
        optimizer_kwargs: dict[str, Any] | None = None,
        store_intermediate_results: bool = False,
        fn_hook: Callable | None = None,
        run_kwargs: dict[str, Any] | None = None,
        filter_kwargs: dict[str, Any] | None = None,
        opt_params: dict[str, Any] | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        additional_fields: dict[str, Any] | None = None,
    ) -> OptSchema:
        """Run a geometry optimization with ASE as the optimizer.

        Parameters
        ----------
        atoms
            Atoms object
        relax_cell
            Whether to relax the cell
        fmax
            Maximum force change
        max_steps
            Maximum number of steps
        optimizer
            Optimizer class to use
        optimizer_kwargs
            Dictionary of kwargs for the optimizer. Takes all valid kwargs for ASE
            Optimizer classes. Refer to `_set_sella_kwargs` for Sella-related
            kwargs and how they are set.
        store_intermediate_results
            Whether to store the files generated at each intermediate step in the
            optimization. If enabled, they will be stored in a directory named
            `stepN` where `N` is the step number, starting at 0.
        fn_hook
            Function to call after each step of the optimization.
        run_kwargs
            Dictionary of kwargs for the `Runner.run_opt` method.
        filter_kwargs
            Dictionary of kwargs for the `FrechetCellFilter` if relax_cell is True.
        opt_params
            DEPRECATED: Dictionary of custom kwargs for the optimization process. For a list
            of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
        copy_files
            Files to copy to runtime directory
        additional_fields
            Additional fields for results

        Returns
        -------
        OptSchema
            Results dictionary
        """
        additional_fields = additional_fields or {}
        additional_fields = {
            "name": f"{self.calculator.__class__.__name__} Relax"
        } | additional_fields

        opt_params = opt_params or {}
        if opt_params:
            warnings.warn(
                "opt_params is deprecated. Use the individual kwargs instead.",
                DeprecationWarning,
                stacklevel=2,
            )

        params = {
            "relax_cell": relax_cell,
            "fmax": fmax,
            "max_steps": max_steps,
            "optimizer": optimizer,
            "optimizer_kwargs": optimizer_kwargs,
            "store_intermediate_results": store_intermediate_results,
            "fn_hook": fn_hook,
            "run_kwargs": run_kwargs,
            "filter_kwargs": filter_kwargs,
        }
        opt_params = params | opt_params

        dyn = Runner(atoms, self.calculator, copy_files=copy_files).run_opt(
            **opt_params
        )
        return Summarize(additional_fields=additional_fields).opt(dyn)

    def vib(
        self,
        atoms: Atoms,
        is_molecule: bool = False,
        vib_kwargs: dict[str, Any] | None = None,
        additional_fields: dict[str, Any] | None = None,
    ) -> VibSchema:
        """Run a vibrational frequency calculation.

        Parameters
        ----------
        atoms
            Atoms object
        vib_kwargs
            Dictionary of kwargs for the [ase.vibrations.Vibrations][] class.
        is_molecule
            Whether the system is a molecular system. If True, the vibrational
            frequencies are cut to 3N-5 or 3N-6.
        additional_fields
            Additional fields for results

        Returns
        -------
        VibSchema
            Results dictionary
        """
        additional_fields = additional_fields or {}
        additional_fields = {
            "name": f"{self.calculator.__class__.__name__} Vibrations"
        } | additional_fields

        vib_kwargs = vib_kwargs or {}
        vib = Runner(atoms, self.calculator).run_vib(vib_kwargs=vib_kwargs)
        return VibSummarize(vib, additional_fields=additional_fields).vib(
            is_molecule=is_molecule
        )
