"""
Materials Project-compatible recipes

This set of recipes is meant to be compatible with the Materials Project
Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np

from quacc import flow, job
from quacc.calculators.vasp import Vasp
from quacc.schemas import fetch_atoms
from quacc.schemas.vasp import summarize_run
from quacc.utils.calc import run_calc
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.vasp import VaspSchema

    class MPRelaxFlowSchema(VaspSchema):
        prerelax: VaspSchema


@job
def mp_prerelax_job(
    atoms: Atoms | dict,
    preset: str | None = "MPScanSet",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to pre-relax a structure with Materials Project settings. By
    default, this uses a PBEsol pre-relax step.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Overrides the following
        defaults:

        ```python
        {"ediffg": -0.05, "xc": "pbesol"}
        ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from quacc.schemas.vasp.summarize_run
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    defaults = {"ediffg": -0.05, "xc": "pbesol"}
    flags = merge_dicts(defaults, calc_swaps, remove_empties=False)

    atoms.calc = Vasp(atoms, preset=preset, **flags)
    atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(atoms, additional_fields={"name": "MP-Prerelax"})


@job
def mp_relax_job(
    atoms: Atoms | dict,
    preset: str | None = "MPScanSet",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to relax a structure with Materials Project settings. By default,
    this uses an r2SCAN relax step.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator. Overrides the following
        defaults: `{}`.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from quacc.schemas.vasp.summarize_run
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    atoms.calc = Vasp(atoms, preset=preset, **calc_swaps)
    atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(atoms, additional_fields={"name": "MP-Relax"})


@flow
def mp_relax_flow(
    atoms: Atoms | dict,
    prerelax: Callable | None = mp_prerelax_job,
    relax: Callable | None = mp_relax_job,
    prerelax_kwargs: dict | None = None,
    relax_kwargs: dict | None = None,
) -> MPRelaxFlowSchema:
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
    MPRelaxFlowSchema
        Dictionary of results
    """
    prerelax_kwargs = prerelax_kwargs or {}
    relax_kwargs = relax_kwargs or {}

    # Run the prerelax
    prerelax_results = prerelax(atoms, **prerelax_kwargs)

    # Update KSPACING arguments
    bandgap = prerelax_results["output"].get("bandgap", 0)
    if bandgap < 1e-4:
        kspacing_swaps = {"kspacing": 0.22, "sigma": 0.2, "ismear": 2}
    else:
        rmin = 25.22 - 2.87 * bandgap
        kspacing = 2 * np.pi * 1.0265 / (rmin - 1.0183)
        kspacing_swaps = {"kspacing": min(kspacing, 0.44), "ismear": -5, "sigma": 0.05}

    relax_kwargs["calc_swaps"] = kspacing_swaps | relax_kwargs.get("calc_swaps", {})

    # Run the relax
    relax_results = relax(prerelax_results, copy_files=["WAVECAR"], **relax_kwargs)
    relax_results["prerelax"] = prerelax_results

    return relax_results
