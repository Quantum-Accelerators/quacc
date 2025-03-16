"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT
from ase.optimize import BFGS

from quacc import job
from quacc.recipes.common.core import Recipe

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from ase.optimize import Dynamics

    from quacc.types import OptSchema, RunSchema


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
    return Recipe(EMT).static(atoms, additional_fields=additional_fields, **calc_kwargs)


@job
def relax_job(
    atoms: Atoms,
    relax_cell: bool = False,
    fmax: float | None = 0.01,
    max_steps: int = 1000,
    optimizer: type[Dynamics] = BFGS,
    optimizer_kwargs: dict[str, Any] | None = None,
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
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        optimizer_kwargs=optimizer_kwargs,
        additional_fields=additional_fields,
        **calc_kwargs,
    )
