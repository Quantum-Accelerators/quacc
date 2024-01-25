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
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.schemas._aliases.vasp import MPRelaxFlowSchema, VaspSchema


@job
def mp_prerelax_job(
    atoms: Atoms,
    preset: str | None = "MPScanSet",
    bandgap: float | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
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
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to `ase.calculators.vasp.vasp.Vasp`.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][].
        See the type-hint for the data structure.
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
    copy_files: str | Path | list[str | Path] | None = None,
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
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to `ase.calculators.vasp.vasp.Vasp`.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][].
        See the type-hint for the data structure.
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
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> MPRelaxFlowSchema:
    """
    Workflow consisting of:

    1. MP-compatible pre-relax
        - name: "mp_prerelax_job"
        - job: [quacc.recipes.vasp.mp.mp_prerelax_job][]

    2. MP-compatible relax
        - name: "mp_relax_job"
        - job: [quacc.recipes.vasp.mp.mp_relax_job][]

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictinoary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    MPRelaxFlowSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    mp_prerelax_job_, mp_relax_job_ = customize_funcs(
        ["mp_prerelax_job", "mp_relax_job"],
        [mp_prerelax_job, mp_relax_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    # Run the prerelax
    prerelax_results = mp_prerelax_job_(atoms)

    # Run the relax
    relax_results = mp_relax_job_(
        prerelax_results["atoms"],
        bandgap=prerelax_results["output"]["bandgap"],
        copy_files=[
            Path(prerelax_results["dir_name"]) / "CHGCAR",
            Path(prerelax_results["dir_name"]) / "WAVECAR",
        ],
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
