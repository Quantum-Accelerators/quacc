"""Molecular dynamics recipes for VASP."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from quacc import job
from quacc.recipes.vasp._base import run_and_summarize

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema
    from quacc.wflow_tools.job_argument import Copy


@job
def md_job(
    atoms: Atoms,
    preset: str | None = "DefaultSetGGA",
    timestep: float = 0.5,
    nsteps: int = 2000,
    ensemble: Literal["nve", "nvt", "npt"] = "nvt",
    temperature: float = 300.0,
    pressure: float = 1.0,
    copy_files: SourceDirectory | Copy | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Run a VASP ab initio molecular dynamics calculation.

    Parameters
    ----------
    atoms
        Atoms object.
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    timestep
        Molecular dynamics timestep in femtoseconds, corresponding to POTIM.
    nsteps
        Number of molecular dynamics steps, corresponding to NSW.
    ensemble
        Thermodynamic ensemble. Supported values are "nve", "nvt", and "npt".
        The value is case-insensitive.
    temperature
        Target temperature in K. For NVE calculations, this temperature is used
        to initialize velocities when velocities are not already present.
    pressure
        External pressure in bar. It will be
        automatically converted to kB. Only used for NPT.
    copy_files
        Files to copy and decompress into the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom keyword arguments for the VASP calculator. These values override
        the recipe defaults. Set a value to None to remove a pre-existing key entirely.
    """
    ensemble = ensemble.lower()

    calc_defaults: dict[str, Any] = {
        "ibrion": 0,
        "isym": 0,
        "nsw": nsteps,
        "potim": timestep,
        "tebeg": temperature,
        "ediffg": None,
    }

    if ensemble == "nvt":
        calc_defaults |= {"isif": 2, "mdalgo": 2, "smass": 0}
    elif ensemble == "npt":
        n_species = len(set(atoms.get_chemical_symbols()))
        calc_defaults |= {
            "isif": 3,
            "mdalgo": 3,
            "langevin_gamma": [10.0] * n_species,
            "langevin_gamma_l": 1.0,
            "pstress": pressure / 1000.0,
        }
    elif ensemble == "nve":
        calc_defaults |= {"isif": 2, "mdalgo": 1, "andersen_prob": 0.0}
    else:
        raise ValueError(
            f"Unsupported ensemble: {ensemble}. Supported ensembles are 'nve', 'nvt', and 'npt'."
        )

    return run_and_summarize(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": f"VASP {ensemble.upper()} Molecular Dynamics"}
        | (additional_fields or {}),
        copy_files=copy_files,
    )
