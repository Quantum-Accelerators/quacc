"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT

from quacc import job
from quacc.recipes._base import BaseRecipe

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import OptParams, OptSchema, RunSchema


class EMTRecipe(BaseRecipe):
    """Base class for EMT recipes."""

    def __init__(self):
        """Initialize EMT recipe.

        Parameters
        ----------
        name
            Optional name override
        """
        super().__init__(EMT)


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
    recipe = EMTRecipe()
    return recipe.run_static(atoms, additional_fields=additional_fields, **calc_kwargs)


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
    recipe = EMTRecipe()
    return recipe.run_relax(
        atoms,
        relax_cell=relax_cell,
        opt_params=opt_params,
        additional_fields=additional_fields,
        **calc_kwargs,
    )
