"""
Materials Project-compatible recipes

This set of recipes is meant to be compatible with the Materials Project
Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801
"""
from __future__ import annotations

import covalent as ct
import numpy as np
from ase import Atoms

from quacc.calculators.vasp import Vasp
from quacc.schemas.atoms import fetch_atoms
from quacc.schemas.vasp import VaspSchema, summarize_run
from quacc.util.calc import run_calc


@ct.electron
def mp_prerelax_job(
    atoms: Atoms | dict,
    preset: str | None = "MPScanSet",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to pre-relax a structure with Materials Project settings.
    By default, this uses a PBEsol pre-relax step.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    preset
        Preset to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from quacc.schemas.vasp.summarize_run
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    defaults = {"ediffg": -0.05, "xc": "pbesol"}
    flags = defaults | calc_swaps

    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(atoms, additional_fields={"name": "MP-Prerelax"})


@ct.electron
def mp_relax_job(
    atoms: Atoms | dict,
    preset: str | None = "MPScanSet",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to relax a structure with Materials Project settings.
    By default, this uses an r2SCAN relax step.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    preset
        Preset to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from quacc.schemas.vasp.summarize_run
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    calc = Vasp(atoms, preset=preset, **calc_swaps)
    atoms.calc = calc
    atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(atoms, additional_fields={"name": "MP-Relax"})


def mp_relax_flow(
    atoms: Atoms | dict,
    prerelax: ct.electron | None = mp_prerelax_job,
    relax: ct.electron | None = mp_relax_job,
    prerelax_kwargs: dict | None = None,
    relax_kwargs: dict | None = None,
) -> VaspSchema:
    """
    Workflow consisting of:

    1. MP-compatible pre-relax

    2. MP-compatible relax

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    prerelax
        Default to use for the pre-relaxation.
    relax
        Default to use for the relaxation.
    prerelax_kwargs
        Additional keyword arguments to pass to the pre-relaxation calculation.
    relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.

    Returns
    -------
    VaspSchema
        Dictionary results from quacc.schemas.vasp.summarize_run
    """
    prerelax_kwargs = prerelax_kwargs or {}
    relax_kwargs = relax_kwargs or {}

    # Run the prerelax
    prerelax_results = prerelax(atoms, **prerelax_kwargs)

    # Update KSPACING arguments
    bandgap = prerelax_results["output"]["bandgap"]
    if bandgap < 1e-4:
        kspacing_swaps = {"kspacing": 0.22, "sigma": 0.2, "ismear": 2, "kpts": None}
    else:
        rmin = 25.22 - 2.87 * bandgap
        kspacing = 2 * np.pi * 1.0265 / (rmin - 1.0183)
        kspacing_swaps = {
            "kspacing": min(kspacing, 0.44),
            "ismear": -5,
            "sigma": 0.05,
            "kpts": None,
        }
    relax_kwargs["calc_swaps"] = kspacing_swaps | relax_kwargs.get("calc_swaps", {})

    # Run the relax
    return relax(prerelax_results, copy_files=["WAVECAR"], **relax_kwargs)
