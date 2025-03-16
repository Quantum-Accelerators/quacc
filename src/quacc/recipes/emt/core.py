"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT

from quacc import job
from quacc.recipes._base import Recipe

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import OptParams, OptSchema, RunSchema


@job
def static_job(
    atoms: Atoms, additional_fields: dict[str, Any] | None = None, **calc_kwargs
) -> RunSchema:
    """Carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    additional_fields
        Additional fields for results
    **calc_kwargs
        Calculator parameters that override defaults

    Returns
    -------
    RunSchema
        Results dictionary
    """
    return Recipe(EMT).static(
        atoms,
        additional_fields={"name": "EMT Static"} | (additional_fields or {}),
        **calc_kwargs,
    )


@job
def relax_job(
    atoms: Atoms,
    relax_cell: bool = False,
    opt_params: OptParams | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    relax_cell
        Whether to relax the cell
    opt_params
        Dictionary of custom kwargs for the optimization process
    additional_fields
        Additional fields for results
    **calc_kwargs
        Calculator parameters that override defaults

    Returns
    -------
    OptSchema
        Results dictionary
    """
    return Recipe(EMT).relax(
        atoms,
        relax_cell=relax_cell,
        opt_params=opt_params,
        additional_fields={"name": "EMT Relax"} | (additional_fields or {}),
        **calc_kwargs,
    )
