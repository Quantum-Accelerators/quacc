"""Core recipes for machine-learned interatomic potentials."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import BFGS

from quacc import job
from quacc.recipes.common.core import Recipe
from quacc.recipes.mlp._defaults import pick_calculator

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms
    from ase.optimize import Dynamics

    from quacc.types import OptSchema, RunSchema


@job
def static_job(
    atoms: Atoms,
    method: Literal["mace-mp-0", "m3gnet", "chgnet", "sevennet", "orb", "fairchem"],
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        Universal ML interatomic potential method to use
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the underlying calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `mace.calculators.mace_mp`, `chgnet.model.dynamics.CHGNetCalculator`,
        `matgl.ext.ase.M3GNetCalculator`, `sevenn.sevennet_calculator.SevenNetCalculator`,
        `orb_models.forcefield.calculator.ORBCalculator`,
        `fairchem.core.common.relaxation.ase_utils.OCPCalculator` calculators.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    calc, calc_defaults, version = pick_calculator(method, **calc_kwargs)
    additional_fields = additional_fields or {}
    additional_fields |= {"mlp_version": version}
    return Recipe(calc, calc_defaults=calc_defaults).static(
        atoms, additional_fields=additional_fields, **calc_kwargs
    )


@job
def relax_job(
    atoms: Atoms,
    method: Literal["mace-mp-0", "m3gnet", "chgnet", "sevennet", "orb", "fairchem"],
    relax_cell: bool = False,
    fmax: float | None = 0.05,
    max_steps: int = 1000,
    optimizer: type[Dynamics] = BFGS,
    optimizer_kwargs: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
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
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the underlying calculator. For a list of available
        keys, refer to the `mace.calculators.mace_mp`, `chgnet.model.dynamics.CHGNetCalculator`,
        `matgl.ext.ase.M3GNetCalculator`, `sevenn.sevennet_calculator.SevenNetCalculator`,
        `orb_models.forcefield.calculator.ORBCalculator`,
        `fairchem.core.common.relaxation.ase_utils.OCPCalculator` calculators.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    calc, calc_defaults, version = pick_calculator(method)
    additional_fields = additional_fields or {}
    additional_fields |= {"mlp_version": version}
    return Recipe(calc, calc_defaults=calc_defaults).relax(
        atoms,
        relax_cell=relax_cell,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        optimizer_kwargs=optimizer_kwargs,
        additional_fields=additional_fields,
        **calc_kwargs,
    )
