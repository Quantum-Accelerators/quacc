"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE

from quacc import job
from quacc.recipes.emt._base import base_job, base_opt_job

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema


@job
def static_job(atoms: Atoms, **calc_kwargs) -> RunSchema:
    """
    Carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.emt.EMT` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """
    calc_defaults = {}

    return base_job(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"EMT Static"},
    )


@job
def relax_job(
    atoms: Atoms,
    relax_cell: bool = False,
    opt_params: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    relax_cell
        Whether to relax the cell
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.emt.EMT` calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in
        [quacc.schemas.ase.summarize_opt_run][]
    """
    calc_defaults = {}
    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}

    return base_opt_job(
        atoms,
        relax_cell=relax_cell,
        opt_defaults=opt_defaults,
        opt_swaps=opt_params,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "EMT Relax"},
    )
