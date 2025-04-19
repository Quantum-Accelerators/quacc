"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT

from quacc import job
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import OptParams, OptSchema, RunSchema


@job
def static_job(
    atoms: Atoms, additional_fields: dict[str, Any] | None = None, **calc_kwargs
) -> RunSchema:
    """
    Carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.emt.EMT][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    calc = EMT(**calc_kwargs)
    final_atoms = Runner(atoms, calc).run_calc()

    return Summarize(
        additional_fields={"name": "EMT Static"} | (additional_fields or {})
    ).run(final_atoms, atoms)


@job
def relax_job(
    atoms: Atoms,
    relax_cell: bool = False,
    opt_params: OptParams | None = None,
    additional_fields: dict[str, Any] | None = None,
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
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.emt.EMT][] calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    opt_params = opt_params or {}

    calc = EMT(**calc_kwargs)
    dyn = Runner(atoms, calc).run_opt(relax_cell=relax_cell, **opt_params)

    return Summarize(
        additional_fields={"name": "EMT Relax"} | (additional_fields or {})
    ).opt(dyn)
