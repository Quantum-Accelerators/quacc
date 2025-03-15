"""Core recipes for DFTB+"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from ase.calculators.dftb import Dftb

from quacc import job
from quacc.recipes._base import BaseRecipe

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, OptParams, RunSchema, SourceDirectory


class DFTBRecipe(BaseRecipe):
    """Base class for DFTB+ recipes."""

    def __init__(self, method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"] = "GFN2-xTB"):
        """Initialize DFTB+ recipe.

        Parameters
        ----------
        name
            Name of the recipe
        method
            Method to use
        """
        calc_defaults = {
            "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
            "Hamiltonian_MaxSccIterations": 200,
        }
        if "xtb" in method.lower():
            calc_defaults["Hamiltonian_Method"] = method

        super().__init__(Dftb, calc_defaults)

    def run_static(
        self,
        atoms: Atoms,
        calc_swaps: dict[str, Any] | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        additional_fields: dict[str, Any] | None = None,
        **calc_kwargs,
    ) -> RunSchema:
        """Run a static calculation.

        Parameters
        ----------
        atoms
            Atoms object
        calc_swaps
            Dictionary of custom kwargs that override defaults
        copy_files
            Files to copy to runtime directory
        additional_fields
            Additional fields for results
        **calc_kwargs
            Additional calculator kwargs

        Returns
        -------
        RunSchema
            Results dictionary
        """
        # Handle k-points
        kpts = calc_kwargs.pop("kpts", None)
        if "kpts" not in self.calc_defaults:
            calc_kwargs["kpts"] = kpts or ((1, 1, 1) if atoms.pbc.any() else None)

        return super().run_static(
            atoms,
            calc_swaps=calc_swaps,
            copy_files=copy_files,
            additional_fields=additional_fields,
            **calc_kwargs,
        )

    def run_relax(
        self,
        atoms: Atoms,
        relax_cell: bool = False,
        opt_params: OptParams | None = None,
        calc_swaps: dict[str, Any] | None = None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        additional_fields: dict[str, Any] | None = None,
        **calc_kwargs,
    ) -> RunSchema:
        """Run a geometry optimization.

        Parameters
        ----------
        atoms
            Atoms object
        relax_cell
            Whether to relax the cell
        opt_params
            Optimization parameters
        calc_swaps
            Dictionary of custom kwargs that override defaults
        copy_files
            Files to copy to runtime directory
        additional_fields
            Additional fields for results
        **calc_kwargs
            Additional calculator kwargs

        Returns
        -------
        RunSchema
            Results dictionary
        """
        # Handle k-points
        kpts = calc_kwargs.pop("kpts", None)
        if "kpts" not in self.calc_defaults:
            calc_kwargs["kpts"] = kpts or ((1, 1, 1) if atoms.pbc.any() else None)

        # Add DFTB-specific relaxation parameters
        calc_kwargs.update(
            {
                "Driver_": "GeometryOptimization",
                "Driver_AppendGeometries": "Yes",
                "Driver_LatticeOpt": "Yes" if relax_cell else "No",
                "Driver_MaxSteps": 2000,
            }
        )

        return super().run_relax(
            atoms,
            relax_cell=relax_cell,
            opt_params=opt_params,
            calc_swaps=calc_swaps,
            copy_files=copy_files,
            additional_fields=additional_fields,
            **calc_kwargs,
        )


@job
def static_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"] = "GFN2-xTB",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    kpts: tuple | list[tuple] | dict | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        Method to use
    copy_files
        Files to copy to runtime directory
    kpts
        k-point grid to use
    additional_fields
        Additional fields for results
    **calc_kwargs
        Additional calculator kwargs

    Returns
    -------
    RunSchema
        Results dictionary
    """
    recipe = DFTBRecipe(method=method)
    return recipe.run_static(
        atoms,
        copy_files=copy_files,
        additional_fields=additional_fields,
        kpts=kpts,
        **calc_kwargs,
    )


@job
def relax_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"] = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict | None = None,
    relax_cell: bool = False,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """Carry out a structure relaxation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        Method to use
    kpts
        k-point grid to use
    relax_cell
        Whether to relax the cell
    copy_files
        Files to copy to runtime directory
    additional_fields
        Additional fields for results
    **calc_kwargs
        Additional calculator kwargs

    Returns
    -------
    RunSchema
        Results dictionary
    """
    recipe = DFTBRecipe(method=method)
    return recipe.run_relax(
        atoms,
        relax_cell=relax_cell,
        copy_files=copy_files,
        additional_fields=additional_fields,
        kpts=kpts,
        **calc_kwargs,
    )
