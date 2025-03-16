"""Base class for all recipes."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

from ase.optimize import BFGS

from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms
    from ase.calculators.calculator import BaseCalculator
    from ase.dynamics import Dynamics

    from quacc.types import Filenames, OptParams, OptSchema, RunSchema, SourceDirectory


class Recipe:
    """Base class for all recipes."""

    def __init__(
        self,
        calculator_class: type[BaseCalculator],
        calc_defaults: dict[str, Any] | None = None,
    ):
        """Initialize the recipe.

        Parameters
        ----------
        calculator_class
            The ASE calculator class to use
        calc_defaults
            Default calculator parameters
        """
        self.calculator_class = calculator_class
        self.calc_defaults = calc_defaults or {}

    def _prepare_calculator(self, **calc_kwargs) -> BaseCalculator:
        """Prepare the calculator with merged parameters.

        Parameters
        ----------
        **calc_kwargs
            Calculator parameters that override defaults

        Returns
        -------
        BaseCalculator
            Configured calculator instance
        """
        calc_flags = recursive_dict_merge(self.calc_defaults, calc_kwargs)
        return self.calculator_class(**calc_flags)

    def calculate(
        self,
        atoms: Atoms,
        geom_file: str | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        additional_fields: dict[str, Any] | None = None,
        **calc_kwargs,
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
        additional_fields |= {"name": f"{self.calculator_class.__name__}"}

        calc = self._prepare_calculator(**calc_kwargs)
        final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc(
            geom_file=geom_file
        )
        return Summarize(additional_fields=additional_fields).run(final_atoms, atoms)

    def static(
        self,
        atoms: Atoms,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        additional_fields: dict[str, Any] | None = None,
        **calc_kwargs,
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
        **calc_kwargs
            Calculator parameters that override defaults

        Returns
        -------
        RunSchema
            Results dictionary
        """
        additional_fields = additional_fields or {}
        additional_fields |= {"name": f"{self.calculator_class.__name__} Static"}

        return self.calculate(
            atoms,
            copy_files=copy_files,
            additional_fields=additional_fields,
            **calc_kwargs,
        )

    def relax(
        self,
        atoms: Atoms,
        relax_cell: bool = False,
        fmax: float | None = 0.01,
        max_steps: int = 1000,
        optimizer: type[Dynamics] = BFGS,
        optimizer_kwargs: dict[str, Any] | None = None,
        store_intermediate_results: bool = False,
        fn_hook: Callable | None = None,
        run_kwargs: dict[str, Any] | None = None,
        filter_kwargs: dict[str, Any] | None = None,
        opt_params: OptParams | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        additional_fields: dict[str, Any] | None = None,
        **calc_kwargs,
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
        **calc_kwargs
            Calculator parameters that override defaults

        Returns
        -------
        OptSchema
            Results dictionary
        """
        additional_fields = additional_fields or {}
        additional_fields |= {"name": f"{self.calculator_class.__name__} Relax"}

        opt_params = opt_params or {}
        if opt_params:
            warnings.warn(
                "opt_params is deprecated. Use the individual kwargs instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            opt_params |= {
                "fmax": fmax,
                "max_steps": max_steps,
                "optimizer": optimizer,
                "optimizer_kwargs": optimizer_kwargs,
                "store_intermediate_results": store_intermediate_results,
                "fn_hook": fn_hook,
                "run_kwargs": run_kwargs,
                "filter_kwargs": filter_kwargs,
            }

        calc = self._prepare_calculator(**calc_kwargs)
        dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(
            relax_cell=relax_cell,
            fmax=fmax,
            max_steps=max_steps,
            optimizer=optimizer,
            optimizer_kwargs=optimizer_kwargs,
            store_intermediate_results=store_intermediate_results,
            fn_hook=fn_hook,
            run_kwargs=run_kwargs,
            filter_kwargs=filter_kwargs,
            **opt_params,
        )
        return Summarize(additional_fields=additional_fields).opt(dyn)
