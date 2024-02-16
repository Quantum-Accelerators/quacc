"""
Molecular Dynamics recipes for EMT.
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

from quacc import job
from quacc.runners.ase import run_md
from quacc.schemas.ase import summarize_md_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import DynSchema


@job
def md_job(
    atoms: Atoms,
    maxwell_boltzmann_params: dict[str, Any] | None = None,
    md_params: dict[str, Any] | None = None,
    **calc_kwargs,
) -> DynSchema:
    """
    Carry out a md calculation.

    Parameters
    ----------
    atoms
        Atoms object
    maxwell_boltzmann_params
        Dictionary of custom kwargs for the initial temperature distribution. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [ase.md.velocitydistribution.MaxwellBoltzmannDistribution][].
        Quacc has two additional parameters `fix_com` and `fix_rot` to fix the center of mass
        and rotation, respectively when setting the initial temperature. The default is
        `{"temperature": None, "fix_com": True, "fix_rot": True}`. If `temperature` is set to
        `None`, the initial temperature will not be set. For an improved chance of
        success when performing a microcanonical ensemble calculation, the initial
        temperature (units of K) should be set to a value twice the target temperature.
    md_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_md][].
    **calc_kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.emt.EMT` calculator.

    Returns
    -------
    DynSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_md_run][].
        See the type-hint for the data structure.
    """

    md_defaults = {"timestep": 1.0, "steps": 500, "dynamics": VelocityVerlet}

    maxwell_boltzmann_defaults = {"temperature": None, "fix_com": True, "fix_rot": True}

    maxwell_boltzmann_params = recursive_dict_merge(
        maxwell_boltzmann_defaults, maxwell_boltzmann_params
    )

    initial_temperature = maxwell_boltzmann_params.pop("temperature", None)
    fix_com = maxwell_boltzmann_params.pop("fix_com", False)
    fix_rot = maxwell_boltzmann_params.pop("fix_rot", False)

    if initial_temperature:
        MaxwellBoltzmannDistribution(
            atoms, temperature_K=initial_temperature, **maxwell_boltzmann_params
        )
        if fix_com:
            Stationary(atoms)
        if fix_rot:
            ZeroRotation(atoms)

    md_flags = recursive_dict_merge(md_defaults, md_params)

    atoms.calc = EMT(**calc_kwargs)

    dyn = run_md(atoms, **md_flags)

    return summarize_md_run(dyn, additional_fields={"name": "EMT Microcanonical"})
