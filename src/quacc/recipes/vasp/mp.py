"""
Materials Project-compatible recipes

This set of recipes is meant to be compatible with the Materials Project
Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from quacc import flow, job
from quacc.recipes.vasp.core import _base_job

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.vasp import VaspSchema

    class MPRelaxFlowSchema(VaspSchema):
        prerelax: VaspSchema


@job
def mp_prerelax_job(
    atoms: Atoms | dict,
    preset: str | None = "MPR2SCANSet",
    bandgap: float = 0.0,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to pre-relax a structure with Materials Project settings. By
    default, this uses a PBEsol pre-relax step.

    ??? Note

        Calculator Defaults:

        ```python
        {"ediffg": -0.05, "xc": "pbesol", "lwave": True, "lcharg": True} | _get_bandgap_swaps(bandgap)
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    bandgap
        Estimate for the bandgap in eV.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    defaults = {
        "ediffg": -0.05,
        "xc": "pbesol",
        "lwave": True,
        "lcharg": True,
    } | _get_bandgap_swaps(bandgap)

    return _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "MP Pre-Relax"},
        copy_files=copy_files,
    )


@job
def mp_relax_job(
    atoms: Atoms | dict,
    preset: str | None = "MPR2SCANSet",
    bandgap: float = 0.0,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to relax a structure with Materials Project settings. By default,
    this uses an r2SCAN relax step.

    ??? Note

        Calculator Defaults:

        ```python
        {"lcharg": True, "lwave": True} | _get_bandgap_swaps(bandgap)
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    bandgap
        Estimate for the bandgap in eV.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    defaults = {"lcharg": True, "lwave": True} | _get_bandgap_swaps(bandgap)
    return _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "MP Relax"},
        copy_files=copy_files,
    )


@job
def mp_static_job(
    atoms: Atoms | dict,
    preset: str | None = "MPR2SCANSet",
    bandgap: float = 0.0,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to run a static calculation with Materials Project settings. By default,
    this uses an r2SCAN static step.

    ??? Note

        Calculator Defaults:

        ```python
        {
        "algo": "fast",
        "nsw": 0,
        "lcharg": True,
        "lwave": False,
        "lreal": False,
        "ismear": -5,
        } | | _get_bandgap_swaps(bandgap)
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    bandgap
        Estimate for the bandgap in eV.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    defaults = {
        "algo": "fast",
        "nsw": 0,
        "lcharg": True,
        "lwave": False,
        "lreal": False,
        "ismear": -5,
    } | _get_bandgap_swaps(bandgap)
    return _base_job(
        atoms,
        preset=preset,
        defaults=defaults,
        calc_swaps=calc_swaps,
        additional_fields={"name": "MP Static"},
        copy_files=copy_files,
    )


@flow
def mp_relax_flow(
    atoms: Atoms | dict,
    prerelax_job_kwargs: dict | None = None,
    relax_job_kwargs: dict | None = None,
) -> MPRelaxFlowSchema:
    """
    Workflow consisting of:

    1. MP-compatible pre-relax

    2. MP-compatible relax

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    prerelax_job_kwargs
        Additional keyword arguments to pass to the pre-relaxation calculation.
    relax_job_kwargs
        Additional keyword arguments to pass to the relaxation calculation.

    Returns
    -------
    MPRelaxFlowSchema
        Dictionary of results
    """
    prerelax_job_kwargs = prerelax_job_kwargs or {}
    relax_job_kwargs = relax_job_kwargs or {}

    # Run the prerelax
    prerelax_results = mp_prerelax_job(atoms, **prerelax_job_kwargs)

    # Update parameters based on bandgap
    bandgap = prerelax_results["output"]["bandgap"]

    # Run the relax
    relax_results = mp_relax_job(
        prerelax_results,
        bandgap=bandgap,
        copy_files=["CHGCAR", "WAVECAR"],
        **relax_job_kwargs,
    )
    relax_results["prerelax"] = prerelax_results

    return relax_results


def _get_bandgap_swaps(bandgap: float) -> dict:
    """
    Get bandgap-related swaps.

    Paramters
    ---------
    bandgap
        The bandgap, in units of eV.

    Returns
    -------
    dict
        Dictionary of swaps.
    """
    if bandgap <= 1e-4:
        kspacing_swaps = {"kspacing": 0.22}
    else:
        rmin = max(1.5, 25.22 - 2.87 * bandgap)
        kspacing = 2 * np.pi * 1.0265 / (rmin - 1.0183)
        kspacing_swaps = {"kspacing": kspacing if 0.22 < kspacing < 0.44 else 0.44}

    if bandgap == 0:
        smearing_swaps = {"ismear": 2, "sigma": 0.2}
    else:
        smearing_swaps = {"ismear": -5, "sigma": 0.05}

    return kspacing_swaps | smearing_swaps
