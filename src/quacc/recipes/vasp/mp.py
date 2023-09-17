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
    preset: str | None = "MPScanSet",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to pre-relax a structure with Materials Project settings. By
    default, this uses a PBEsol pre-relax step.

    ??? Note

        Calculator Defaults:

        ```python
        {"ediffg": -0.05, "xc": "pbesol"}
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    defaults = {"ediffg": -0.05, "xc": "pbesol"}
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
    preset: str | None = "MPScanSet",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> VaspSchema:
    """
    Function to relax a structure with Materials Project settings. By default,
    this uses an r2SCAN relax step.

    ??? Note

        Calculator Defaults:

        ```python
        {}
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    preset
        Preset to use.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    return _base_job(
        atoms,
        preset=preset,
        defaults={},
        calc_swaps=calc_swaps,
        additional_fields={"name": "MP Relax"},
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

    # Update KSPACING arguments
    bandgap = prerelax_results["output"].get("bandgap", 0)
    if bandgap < 1e-4:
        kspacing_swaps = {"kspacing": 0.22, "sigma": 0.2, "ismear": 2}
    else:
        rmin = 25.22 - 2.87 * bandgap
        kspacing = 2 * np.pi * 1.0265 / (rmin - 1.0183)
        kspacing_swaps = {"kspacing": min(kspacing, 0.44), "ismear": -5, "sigma": 0.05}

    relax_job_kwargs["calc_swaps"] = kspacing_swaps | relax_job_kwargs.get(
        "calc_swaps", {}
    )

    # Run the relax
    relax_results = mp_relax_job(
        prerelax_results, copy_files=["WAVECAR"], **relax_job_kwargs
    )
    relax_results["prerelax"] = prerelax_results

    return relax_results
