"""
Materials Project-compatible recipes.

This set of recipes is meant to be compatible with the Materials Project

!!! Info

    The one true source of Materials Project workflows is
    [atomate2](https://github.com/materialsproject/atomate2).
    If you need an MP-compatible workflow, we strongly encourage you to
    use atomate2 to ensure that all of your settings are fully compatible
    and up-to-date. This module is a best effort to be used at your own
    discretion.
"""

from __future__ import annotations

from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING

from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet, MPStaticSet

from quacc import flow, job
from quacc.recipes.vasp._base import base_fn
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.schemas._aliases.vasp import (
        DoubleRelaxSchema,
        MPGGARelaxFlowSchema,
        MPMetaGGARelaxFlowSchema,
        VaspSchema,
    )


@job
def mp_gga_relax_job(
    atoms: Atoms, copy_files: str | Path | list[str | Path] | None = None, **calc_kwargs
) -> DoubleRelaxSchema:
    """
    Function to (double) relax a structure with the original Materials Project GGA(+U) settings.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [ase.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    DoubleRelaxSchema
        Dictionary of results.
    """

    def _relax(
        atoms: Atoms,
        copy_files: str | Path | list[str | Path] | None = None,
        calc_kwargs: dict[str, Any] | None = None,
    ) -> VaspSchema:
        """A helper function to run a relaxation with the MP GGA settings."""
        calc_defaults = {"pmg_input_set": MPRelaxSet}
        return base_fn(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            additional_fields={"name": "MP GGA Relax"},
            copy_files=copy_files,
        )

    summary1 = _relax(atoms, copy_files=copy_files, **calc_kwargs)
    summary2 = _relax(
        atoms,
        copy_files=[
            Path(summary1["dir_name"]) / "CHGCAR",
            Path(summary1["dir_name"]) / "WAVECAR",
        ],
        **calc_kwargs,
    )

    return {"relax1": summary1, "relax2": summary2}


@job
def mp_gga_static_job(
    atoms: Atoms,
    bandgap: float | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Function to run a static calculation on a structure with the original Materials Project GGA(+U) settings.

    Parameters
    ----------
    atoms
        Atoms object
    bandgap
        The bandgap in eV, if known from a prior calculation.
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [ase.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][].
    """

    calc_defaults = {
        "pmg_input_set": partial(
            MPStaticSet, bandgap=bandgap, small_gap_multiply=[1e-4, 3.125]
        ),
        "algo": "fast",
        "lwave": True,  # Deviation from MP (but logical)
        "lreal": False,
    }
    return base_fn(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP GGA Static"},
        copy_files=copy_files,
    )


@job
def mp_metagga_prerelax_job(
    atoms: Atoms,
    bandgap: float | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Function to pre-relax a structure with Materials Project r2SCAN workflow settings. By default, this
    uses a PBEsol pre-relax step.

    Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801

    Parameters
    ----------
    atoms
        Atoms object
    bandgap
        Estimate for the bandgap in eV.
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [ase.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.vasp_summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "pmg_input_set": partial(
            MPScanRelaxSet, bandgap=bandgap or 0.0, auto_ismear=False
        ),
        "ediffg": -0.05,
        "gga": "PS",
        "laechg": False,  # Deviation from MP (but logical)
        "lvtot": False,  # Deviation from MP (but logical)
        "lwave": True,
        "metagga": None,
    }
    return base_fn(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP Meta-GGA Pre-Relax"},
        copy_files=copy_files,
    )


@job
def mp_metagga_relax_job(
    atoms: Atoms,
    bandgap: float | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
    **calc_kwargs,
) -> DoubleRelaxSchema:
    """
    Function to (double) relax a structure with Materials Project r2SCAN workflow settings. By default, this uses
    an r2SCAN relax step.

    Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801

    Parameters
    ----------
    atoms
        Atoms object
    bandgap
        Estimate for the bandgap in eV.
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [ase.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    DoubleRelaxSchema
        Dictionary of results.
    """

    def _relax(
        atoms: Atoms,
        copy_files: str | Path | list[str | Path] | None = None,
        bandgap: float | None = None,
        calc_kwargs: dict[str, Any] | None = None,
    ) -> VaspSchema:
        """A helper function to run a relaxation with the MP r2SCAN settings."""
        calc_defaults = {
            "pmg_input_set": partial(
                MPScanRelaxSet, bandgap=bandgap or 0.0, auto_ismear=False
            ),
            "laechg": False,  # Deviation from MP (but logical)
            "lvtot": False,  # Deviation from MP (but logical)
            "lwave": True,
        }
        return base_fn(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            additional_fields={"name": "MP Meta-GGA Relax"},
            copy_files=copy_files,
        )

    summary1 = _relax(atoms, copy_files=copy_files, bandgap=bandgap, **calc_kwargs)
    summary2 = _relax(
        atoms,
        copy_files=[
            Path(summary1["dir_name"]) / "CHGCAR",
            Path(summary1["dir_name"]) / "WAVECAR",
        ],
        bandgap=bandgap,
        **calc_kwargs,
    )

    return {"relax1": summary1, "relax2": summary2}


@job
def mp_metagga_static_job(
    atoms: Atoms,
    bandgap: float | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Function to run a static calculation on a structure with r2SCAN workflow Materials Project settings.
    By default, this uses an r2SCAN static step.

    Parameters
    ----------
    atoms
        Atoms object
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

    calc_defaults = {
        "pmg_input_set": partial(
            MPScanRelaxSet, bandgap=bandgap or 0.0, auto_ismear=False
        ),
        "algo": "fast",
        "ismear": -5,
        "lreal": False,
        "lwave": True,  # Deviation from MP (but logical)
        "nsw": 0,
    }
    return base_fn(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP Meta-GGA Static"},
        copy_files=copy_files,
    )


@flow
def mp_gga_relax_flow(
    atoms: Atoms,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> MPGGARelaxFlowSchema:
    """
    Materials Project GGA workflow consisting of:

    1. MP-compatible (double) relax
        - name: "mp_gga_relax_job"
        - job: [quacc.recipes.vasp.mp.mp_gga_relax_job][]

    2. MP-compatible static
        - name: "mp_gga_static_job"
        - job: [quacc.recipes.vasp.mp.mp_gga_static_job][]

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
    MPGGARelaxFlowSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    (mp_gga_relax_job_, mp_gga_static_job_) = customize_funcs(
        ["mp_gga_relax_job", "mp_gga_static_job"],
        [mp_gga_relax_job, mp_gga_static_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    # Run the relax
    relax_results = mp_gga_relax_job_(atoms)

    # Run the static
    static_results = mp_gga_static_job_(
        relax_results["relax2"]["atoms"],
        bandgap=relax_results["relax2"]["output"]["bandgap"],
        copy_files=[
            Path(relax_results["relax2"]["dir_name"]) / "CHGCAR",
            Path(relax_results["relax2"]["dir_name"]) / "WAVECAR",
        ],
    )

    return {"relax": relax_results, "static": static_results}


@flow
def mp_metagga_relax_flow(
    atoms: Atoms,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> MPMetaGGARelaxFlowSchema:
    """
    Materials Project r2SCAN workflow consisting of:

    1. MP-compatible pre-relax
        - name: "mp_metagga_prerelax_job"
        - job: [quacc.recipes.vasp.mp.mp_metagga_prerelax_job][]

    2. MP-compatible (double) relax
        - name: "mp_metagga_relax_job"
        - job: [quacc.recipes.vasp.mp.mp_metagga_relax_job][]

    3. MP-compatible static
        - name: "mp_metagga_static_job"
        - job: [quacc.recipes.vasp.mp.mp_metagga_static_job][]

    Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801

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
    MPMetaGGARelaxFlowSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    (
        mp_metagga_prerelax_job_,
        mp_metagga_relax_job_,
        mp_metagga_static_job_,
    ) = customize_funcs(
        ["mp_metagga_prerelax_job", "mp_metagga_relax_job", "mp_metagga_static_job"],
        [mp_metagga_prerelax_job, mp_metagga_relax_job, mp_metagga_static_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    # Run the prerelax
    prerelax_results = mp_metagga_prerelax_job_(atoms)

    # Run the relax
    relax_results = mp_metagga_relax_job_(
        prerelax_results["atoms"],
        bandgap=prerelax_results["output"]["bandgap"],
        copy_files=[
            Path(prerelax_results["dir_name"]) / "CHGCAR",
            Path(prerelax_results["dir_name"]) / "WAVECAR",
        ],
    )

    # Run the static
    static_results = mp_metagga_static_job_(
        relax_results["relax2"]["atoms"],
        bandgap=relax_results["relax2"]["output"]["bandgap"],
        copy_files=[
            Path(relax_results["relax2"]["dir_name"]) / "CHGCAR",
            Path(relax_results["relax2"]["dir_name"]) / "WAVECAR",
        ],
    )

    return {
        "prerelax": prerelax_results,
        "relax": relax_results,
        "static": static_results,
    }
