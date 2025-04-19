"""Core recipes for DFTB+"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from ase.calculators.dftb import Dftb

from quacc import job
from quacc.recipes.common.core import Recipe
from quacc.recipes.dftb._defaults import _GEOM_FILE, create_dftb_defaults
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, RunSchema, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"] = "GFN2-xTB",
    kpts: tuple | list[tuple] | dict | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
    kpts
        k-point grid to use
    copy_files
        Files to copy to runtime directory
    additional_fields
        Metadata to store in the results
    **calc_kwargs
        Calculator parameters to pass to [ase.calculators.dftb.Dftb][]

    Returns
    -------
    RunSchema
        Results dictionary
    """
    calc_defaults = create_dftb_defaults(
        method=method, kpts=kpts, is_periodic=bool(atoms.pbc.any())
    )
    calc_params = recursive_dict_merge(calc_defaults, calc_kwargs)
    calc = Dftb(**calc_params)
    return Recipe(calc).run(
        atoms,
        geom_file=_GEOM_FILE,
        copy_files=copy_files,
        additional_fields=additional_fields,
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
        Metadata to store in the results
    **calc_kwargs
        Calculator parameters to pass to [ase.calculators.dftb.Dftb][]

    Returns
    -------
    RunSchema
        Results dictionary
    """
    calc_defaults = create_dftb_defaults(
        method=method, kpts=kpts, is_periodic=bool(atoms.pbc.any())
    )
    calc_defaults |= {
        "Driver_": "GeometryOptimization",
        "Driver_AppendGeometries": "Yes",
        "Driver_LatticeOpt": "Yes" if relax_cell else "No",
        "Driver_MaxSteps": 2000,
    }
    calc_params = recursive_dict_merge(calc_defaults, calc_kwargs)
    calc = Dftb(**calc_params)
    return Recipe(calc).run(
        atoms,
        geom_file=_GEOM_FILE,
        copy_files=copy_files,
        additional_fields=additional_fields,
    )
