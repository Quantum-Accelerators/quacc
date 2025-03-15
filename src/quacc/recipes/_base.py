"""Base class for all recipes."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from ase.calculators.calculator import BaseCalculator

    from quacc.types import Filenames, OptParams, OptSchema, RunSchema, SourceDirectory


class BaseRecipe:
    """Base class for all recipes."""

    def __init__(
        self,
        calculator_class: type[BaseCalculator],
        name: str | None = None,
        calc_defaults: dict[str, Any] | None = None,
    ):
        """Initialize the recipe.

        Parameters
        ----------
        calculator_class
            The ASE calculator class to use
        name
            Name of the recipe (e.g. "DFTB+ Static"). If not provided,
            will be derived from the calculator class name.
        calc_defaults
            Default calculator parameters
        """
        self.calculator_class = calculator_class
        self.name = name or calculator_class.__name__
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

    def run_static(
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
        calc = self._prepare_calculator(**calc_kwargs)
        final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()
        return Summarize(
            additional_fields={"name": f"{self.name} Static"}
            | (additional_fields or {})
        ).run(final_atoms, atoms)

    def run_relax(
        self,
        atoms: Atoms,
        relax_cell: bool = False,
        opt_params: OptParams | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        additional_fields: dict[str, Any] | None = None,
        **calc_kwargs,
    ) -> OptSchema:
        """Run a geometry optimization.

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
        calc = self._prepare_calculator(**calc_kwargs)
        opt_params = opt_params or {}
        opt_params.setdefault("relax_cell", relax_cell)
        dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(**opt_params)
        return Summarize(
            additional_fields={"name": f"{self.name} Relax"} | (additional_fields or {})
        ).opt(dyn)
