"""
Materials Project-compatible recipes.

This set of recipes is meant to be compatible with the Materials Project
Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801

!!! Info

    The one true source of Materials Project workflows is
    [atomate2](https://github.com/materialsproject/atomate2).
    If you need an MP-compatible workflow, we strongly encourage you to
    use atomate2 to ensure that all of your settings are fully compatible
    and up-to-date. This module is a best effort to be used at your own
    discretion.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from quacc import flow, job
from quacc.recipes.vasp._base import base_fn

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.vasp import MPRelaxFlowSchema, VaspSchema


@job
def mp_prerelax_job(
    atoms: Atoms,
    preset: str | None = "MPScanSet",
    bandgap: float | None = None,
    copy_files: list[str] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Function to pre-relax a structure with Materials Project settings. By default, this
    uses a PBEsol pre-relax step.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    bandgap
        Estimate for the bandgap in eV.
    copy_files
        Files to copy to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `quacc.calculators.vasp.vasp.Vasp` calculator.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    calc_defaults = {
        "ediffg": -0.05,
        "xc": "pbesol",
        "lwave": True,
        "lcharg": True,
    } | _get_bandgap_swaps(bandgap)

    return base_fn(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP Pre-Relax"},
        copy_files=copy_files,
    )


@job
def mp_relax_job(
    atoms: Atoms,
    preset: str | None = "MPScanSet",
    bandgap: float | None = None,
    copy_files: list[str] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Function to relax a structure with Materials Project settings. By default, this uses
    an r2SCAN relax step.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    bandgap
        Estimate for the bandgap in eV.
    copy_files
        Files to copy to the runtime directory.
    **calc_kwargs
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `quacc.calculators.vasp.vasp.Vasp` calculator.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    calc_defaults = {"lcharg": True, "lwave": True} | _get_bandgap_swaps(bandgap)
    return base_fn(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP Relax"},
        copy_files=copy_files,
    )


@flow
def mp_relax_flow(
    atoms: Atoms,
    prerelax_job_kwargs: dict[str, Any] | None = None,
    relax_job_kwargs: dict[str, Any] | None = None,
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
        Additional keyword arguments to pass to [quacc.recipes.vasp.mp.mp_prerelax_job][].
    relax_job_kwargs
        Additional keyword arguments to pass to [quacc.recipes.vasp.mp.mp_relax_job][].

    Returns
    -------
    MPRelaxFlowSchema
        Dictionary of results
    """
    prerelax_job_kwargs = prerelax_job_kwargs or {}
    relax_job_kwargs = relax_job_kwargs or {}

    # Run the prerelax
    prerelax_results = mp_prerelax_job(atoms, **prerelax_job_kwargs)

    # Run the relax
    relax_results = mp_relax_job(
        prerelax_results["atoms"],
        bandgap=prerelax_results["output"]["bandgap"],
        copy_files=[
            Path(prerelax_results["dir_name"]) / "CHGCAR",
            Path(prerelax_results["dir_name"]) / "WAVECAR",
        ],
        **relax_job_kwargs,
    )
    relax_results["prerelax"] = prerelax_results

    return relax_results


def _get_bandgap_swaps(bandgap: float | None = None) -> dict[str, float]:
    """
    Get bandgap-related swaps.

    Parameters
    ---------
    bandgap
        The bandgap, in units of eV.

    Returns
    -------
    dict
        Dictionary of swaps.
    """

    if bandgap is None:
        return {"kspacing": 0.22, "ismear": 0, "sigma": 0.05}
    if bandgap <= 1e-4:
        return {"kspacing": 0.22, "ismear": 2, "sigma": 0.2}
    rmin = max(1.5, 25.22 - 2.87 * bandgap)
    kspacing = 2 * np.pi * 1.0265 / (rmin - 1.0183)
    return {"kspacing": min(kspacing, 0.44), "ismear": -5, "sigma": 0.05}
