"""
Molecular Dynamics recipes for EMT.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT
from ase.md.verlet import VelocityVerlet

from quacc import job
from quacc.runners.ase import Runner
from quacc.schemas.ase import summarize_md_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from ase.md.md import MolecularDynamics

    from quacc.runners.ase import MDParams
    from quacc.schemas._aliases.ase import DynSchema
    from quacc.utils.files import Filenames, SourceDirectory


@job
def md_job(
    atoms: Atoms,
    dynamics: MolecularDynamics = VelocityVerlet,
    dynamics_kwargs: dict[str, Any] | None = None,
    initial_temperature: float | None = None,
    md_params: MDParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> DynSchema:
    """
    Carry out a Molecular Dynamics calculation.

    !!! Important "Units"

        Quacc does not follow ASE standards for Molecular Dynamics units and instead
        uses the following units:

        - Time: femtoseconds (fs)
        - Pressure: GPa
        - Temperature: Kelvin (K)
        - Compressibility: 1/GPa

    Parameters
    ----------
    atoms
        Atoms object
    dynamics
        ASE `MolecularDynamics` class to use, from `ase.md.md.MolecularDynamics`.
    dynamics_kwargs
        Dictionary of custom keyword arguments for the chosen dynamics class.
    initial_temperature
        Initial temperature (in K) to specify via a Maxwell-Boltzmann distribution.
    md_params
        Dictionary of custom kwargs for the MD run. For a list of available
        keys, refer to [quacc.runners.ase.Runner.run_md][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
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
    md_defaults = {
        "dynamics": dynamics,
        "dynamics_kwargs": dynamics_kwargs,
        "maxwell_boltzmann_kwargs": {"temperature_K": initial_temperature}
        if initial_temperature
        else None,
    }
    md_params = recursive_dict_merge(md_defaults, md_params)

    calc = EMT(**calc_kwargs)
    dyn = Runner(atoms, calc, copy_files=copy_files).run_md(**md_params)

    return summarize_md_run(dyn, additional_fields={"name": "EMT MD"})
