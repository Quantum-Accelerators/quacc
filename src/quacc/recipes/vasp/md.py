from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.calculators.vasp import Vasp
from quacc.recipes.vasp._base import run_and_summarize

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from quacc.types import VaspSchema

@job
def md_job(
    atoms: Atoms,
    timestep: float = 1.0,
    ensemble: str = "NVT",
    temperature: float = 300.0,
    pressure: float = 0.0,
    copy_files: list[str] | None = None,
    **kwargs,
) -> VaspSchema:
    """
    VASP Ab Initio Molecular Dynamics (AIMD) job.

    Parameters
    ----------
    atoms
        Atoms object
    timestep
        Timestep in femtoseconds (POTIM).
    ensemble
        Thermodynamic ensemble: "NVT", "NPT", or "NVE".
    temperature
        Target temperature in Kelvin (TEBEG).
    pressure
        Target pressure in kB (PSTRESS). Only used if ensemble is NPT.
    copy_files
        Files to copy to the runtime directory.
    **kwargs
        Custom kwargs for the Vasp calculator.

    Returns
    -------
    VaspSchema
        Dictionary of results from quacc.schemas.vasp.vasp_summarize_run
    """
    # Core AIMD parameters
    user_incar_settings = {
        "IBRION": 0,
        "POTIM": timestep,
        "TEBEG": temperature,
    }

    # Set ensemble-specific parameters
    ens = ensemble.upper()
    if ens == "NVT":
        user_incar_settings.update({"MDALGO": 2, "ISIF": 2})
    elif ens == "NPT":
        user_incar_settings.update({"MDALGO": 3, "ISIF": 3, "PSTRESS": pressure})
    elif ens == "NVE":
        user_incar_settings.update({"MDALGO": 0, "ISIF": 2})
    else:
        raise ValueError(f"Unknown ensemble: {ensemble}. Choose from NVT, NPT, or NVE.")

    # Merge custom user kwargs
    user_incar_settings.update(kwargs)

    # Setup VASP calculator
    calc = Vasp(atoms, preset="AIMDSet", incar_copier=False, **user_incar_settings)
    atoms.calc = calc

    # Run calculation and return summary
    return run_and_summarize(atoms, copy_files=copy_files)

