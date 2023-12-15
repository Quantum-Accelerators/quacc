"""Core recipes for universal machine-learned interatomic potentials."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE

from quacc import job
from quacc.recipes.mlp._base import _pick_calculator
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema


@job
def static_job(
    atoms: Atoms, method: Literal["umace", "m3gnet", "chgnet"], **calc_kwargs
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        Universal ML interatomic potential method to use
    **calc_kwargs
        Custom kwargs for the underlying calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `mace.calculators.mace_mp`, `chgnet.model.dynamics.CHGNetCalculator`,
        or `matgl.ext.ase.M3GNetCalculator` calculators.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    atoms.calc = _pick_calculator(method, **calc_kwargs)
    final_atoms = run_calc(atoms)
    return summarize_run(
        final_atoms, input_atoms=atoms, additional_fields={"name": f"{method} Static"}
    )


@job
def relax_job(
    atoms: Atoms,
    method: Literal["umace", "m3gnet", "chgnet"],
    relax_cell: bool = False,
    opt_params: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object
    method
        Universal ML interatomic potential method to use
    relax_cell
        Whether to relax the cell.
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    **calc_kwargs
        Custom kwargs for the underlying calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `mace.calculators.mace_mp`, `chgnet.model.dynamics.CHGNetCalculator`,
        or `matgl.ext.ase.M3GNetCalculator` calculators.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """

    atoms.calc = _pick_calculator(method, **calc_kwargs)

    opt_defaults = {"fmax": 0.1, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = merge_dicts(opt_defaults, opt_params)

    dyn = run_opt(atoms, relax_cell=relax_cell, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": f"{method} Relax"})
