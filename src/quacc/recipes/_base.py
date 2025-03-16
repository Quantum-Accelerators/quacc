"""Base class for all recipes."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize, VibSummarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from ase.calculators.calculator import BaseCalculator

    from quacc.types import (
        Filenames,
        OptParams,
        OptSchema,
        RunSchema,
        SourceDirectory,
        VibKwargs,
        VibThermoSchema,
    )


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
        additional_fields.setdefault("name", f"{self.calculator_class.__name__}")

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
        additional_fields.setdefault("name", f"{self.calculator_class.__name__} Static")

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
        opt_params
            Optimization parameters
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
        additional_fields.setdefault("name", f"{self.calculator_class.__name__} Relax")

        opt_params = opt_params or {}
        opt_params.setdefault("relax_cell", relax_cell)

        calc = self._prepare_calculator(**calc_kwargs)
        dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(**opt_params)
        return Summarize(additional_fields=additional_fields).opt(dyn)
