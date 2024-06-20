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

from quacc import job
from quacc.runners.ase import Runner
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
    Carry out a Molecular Dynamics calculation.

    !!! Note "Units"

        Quacc does not follow ASE standards for Molecular Dynamics units.

        Quacc ALWAYS uses the following units:

        - Time: femtoseconds (fs)
        - Pressure: GPa
        - Temperature: Kelvin (K)
        - Compressibility: 1/GPa

    !!! Note "Keywords"

        Additionally, Quacc uses the keywords `fix_com` and `fix_rot` instead of
        `fixcm` and `fixrot` to fix the center of mass and rotation, respectively.

        Shared keywords between various dynamics type such as `timestep` and `steps` should
        be specified in the `md_params` dictionary. Keywords specific to the dynamics type
        should be specified in a dictionary `dynamics_kwargs` inside `md_params`.

        The dynamics type can be specified in `md_params` as `dynamics`.

    Parameters
    ----------
    atoms
        Atoms object
    maxwell_boltzmann_params
        Dictionary of custom kwargs for the initial temperature distribution. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [ase.md.velocitydistribution.MaxwellBoltzmannDistribution][].
        Quacc has two additional parameters `fix_com` and `fix_rot` to fix the center of mass
        and rotation from the initial velocity distribution. Defaults are
        `{"temperature": None, "fix_com": True, "fix_rot": True}`. If `temperature` is set to
        `None`, the initial temperature will not be set.
    md_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.Runner.run_md][].
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

    maxwell_boltzmann_defaults = {"fix_com": True, "fix_rot": True, "temperature": None}
    md_params = md_params or {}

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

    calc = EMT(**calc_kwargs)
    dyn = Runner(atoms, calc).run_md(**md_params)

    return summarize_md_run(dyn, additional_fields={"name": "EMT MD"})
