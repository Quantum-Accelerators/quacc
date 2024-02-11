"""
Core recipes for EMT.

NOTE: This set of minimal recipes is mainly for demonstration purposes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation,
)
from ase.md.verlet import VelocityVerlet
from ase.optimize import FIRE
from ase.units import fs

from quacc import job
from quacc.runners.ase import run_calc, run_md, run_opt
from quacc.schemas.ase import summarize_opt_run, summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import DynSchema, OptSchema, RunSchema


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
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.emt.EMT` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """
    atoms.calc = EMT(**calc_kwargs)
    final_atoms = run_calc(atoms)

    return summarize_run(final_atoms, atoms, additional_fields={"name": "EMT Static"})


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
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.emt.EMT` calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_opt_run][].
        See the type-hint for the data structure.
    """
    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    atoms.calc = EMT(**calc_kwargs)

    dyn = run_opt(atoms, relax_cell=relax_cell, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "EMT Relax"})


@job
def microcanonical_job(
    atoms: Atoms,
    initial_temperature_params: float = None,
    md_params: dict[str, Any] | None = None,
    **calc_kwargs,
) -> DynSchema:
    """
    Carry out a microcanonical ensemble calculation.

    Parameters
    ----------
    atoms
        Atoms object
    relax_cell
        Whether to relax the cell
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.emt.EMT` calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_opt_run][].
        See the type-hint for the data structure.
    """

    md_defaults = {"timestep": 1.0 * fs, "max_steps": 500, "dynamics": VelocityVerlet}

    initial_temperature_defaults = {"temperature": None, "fixcm": True, "fixrot": True}

    initial_temperature_params = recursive_dict_merge(
        initial_temperature_defaults, initial_temperature_params
    )

    temperature = initial_temperature_params.get("temperature")

    if temperature:
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
        if initial_temperature_params.get("fixcm"):
            Stationary(atoms)
        if initial_temperature_params.get("fixrot"):
            ZeroRotation(atoms)

    md_flags = recursive_dict_merge(md_defaults, md_params)

    atoms.calc = EMT(**calc_kwargs)

    dyn = run_md(atoms, **md_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "EMT Relax"})
