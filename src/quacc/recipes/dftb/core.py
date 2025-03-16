"""Core recipes for DFTB+"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from ase.calculators.dftb import Dftb

from quacc import job
from quacc.recipes._base import Recipe

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, RunSchema, SourceDirectory

_GEOM_FILE = "geo_end.gen"


def _create_dftb_defaults(
    method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"],
    kpts: tuple | list[tuple] | dict | None = None,
    is_periodic: bool = True,
) -> dict[str, Any]:
    """Create the default calculator kwargs for DFTB+.

    Parameters
    ----------
    method
        Method to use
    kpts
        k-point grid to use
    is_periodic
        Whether the system is periodic

    Returns
    -------
    dict[str, Any]
        Default calculator kwargs
    """
    calc_defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_MaxSccIterations": 200,
        "kpts": kpts or ((1, 1, 1) if is_periodic else None),
    }
    if "xtb" in method.lower():
        calc_defaults["Hamiltonian_Method"] = method
    return calc_defaults


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
    calc_defaults = _create_dftb_defaults(
        method, kpts=kpts, is_periodic=atoms.pbc.any()
    )
    recipe = Recipe(Dftb, calc_defaults)
    return recipe.calculate(
        atoms,
        geom_file=_GEOM_FILE,
        copy_files=copy_files,
        additional_fields=additional_fields,
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
    calc_defaults = _create_dftb_defaults(
        method, kpts=kpts, is_periodic=atoms.pbc.any()
    )
    calc_defaults.update(
        {
            "Driver_": "GeometryOptimization",
            "Driver_AppendGeometries": "Yes",
            "Driver_LatticeOpt": "Yes" if relax_cell else "No",
            "Driver_MaxSteps": 2000,
        }
    )
    recipe = Recipe(Dftb, calc_defaults)
    return recipe.calculate(
        atoms,
        geom_file=_GEOM_FILE,
        copy_files=copy_files,
        additional_fields=additional_fields,
        **calc_kwargs,
    )
