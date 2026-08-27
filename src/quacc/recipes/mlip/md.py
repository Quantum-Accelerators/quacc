"""Molecular dynamics recipes for universal machine-learned interatomic potentials."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.md.verlet import VelocityVerlet
from ase.units import bar, fs

from quacc import Remove, job
from quacc.recipes.mlip._base import pick_calculator
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms
    from ase.md.md import MolecularDynamics

    from quacc.types import DynSchema, MDParams


@job
def md_job(
    atoms: Atoms,
    library: Literal["fairchem", "matcalc", "rootstock"],
    dynamics: MolecularDynamics = VelocityVerlet,
    steps: int = 1000,
    timestep_fs: float = 1.0,
    temperature_K: float | None = None,
    pressure_bar: float | None = None,
    md_params: MDParams | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> DynSchema:
    """
    Carry out a Molecular Dynamics calculation.

    Parameters
    ----------
    atoms
        Atoms object
    library
        MLIP library to use:
        - `fairchem` passes `**calc_kwargs` to `FAIRChemCalculator.from_model_checkpoint()`
        - `matcalc` passes `**calc_kwargs` to `matcalc.load_fp()`
        - `rootstock` passes `**calc_kwargs` to `rootstock.RootstockCalculator()`
    dynamics
        ASE `MolecularDynamics` class to use, from `ase.md.md.MolecularDynamics`.
    steps
        Number of MD steps to run.
    timestep_fs
        Time step in fs.
    temperature_K
        Temperature in K, if applicable for the given ensemble.
    pressure_bar
        Pressure in bar, if applicable for the given ensemble.
    md_params
        Dictionary of custom kwargs for the MD run. For a list of available
        keys, refer to [quacc.runners.ase.Runner.run_md][].
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the underlying MLIP library.

    Returns
    -------
    DynSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.md][].
        See the type-hint for the data structure.
    """
    md_defaults = {
        "steps": steps,
        "dynamics_kwargs": {
            "timestep": timestep_fs * fs,
            "temperature_K": temperature_K if temperature_K else Remove,
            "pressure_au": pressure_bar * bar if pressure_bar else Remove,
        },
    }
    md_params = recursive_dict_merge(md_defaults, md_params)

    calc = pick_calculator(library, **calc_kwargs)
    dyn = Runner(atoms, calc).run_md(dynamics, **md_params)

    return Summarize(
        additional_fields={"name": f"{library} MD"} | (additional_fields or {})
    ).md(dyn)
