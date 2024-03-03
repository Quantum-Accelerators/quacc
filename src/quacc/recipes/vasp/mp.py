"""
Materials Project-compatible recipes.

This set of recipes is meant to be compatible with the Materials Project

!!! Important

    Make sure that you use the Materials Project-compatible pseudpotential
    versions. The GGA workflows use the old (no version) PAW PBE potentials.
    The meta-GGA workflows currently use the v.54 PAW PBE potentials.

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
from typing import TYPE_CHECKING

from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet, MPStaticSet

from quacc import flow, job
from quacc.recipes.vasp._base import base_fn
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.schemas._aliases.vasp import (
        MPGGARelaxFlowSchema,
        MPMetaGGARelaxFlowSchema,
        VaspSchema,
    )
    from quacc.utils.files import Filenames, SourceDirectory


@job
def mp_gga_relax_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Function to relax a structure with the original Materials Project GGA(+U) settings.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VaspSchema
        Dictionary of results.
    """

    calc_defaults = {"pmg_input_set": MPRelaxSet}
    return base_fn(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        report_mp_corrections=True,
        additional_fields={"name": "MP GGA Relax"},
        copy_files=copy_files,
    )


@job
def mp_gga_static_job(
    atoms: Atoms,
    bandgap: float | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

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
        report_mp_corrections=True,
        additional_fields={"name": "MP GGA Static"},
        copy_files=copy_files,
    )


@job
def mp_metagga_prerelax_job(
    atoms: Atoms,
    bandgap: float | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

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
        report_mp_corrections=True,
        additional_fields={"name": "MP Meta-GGA Pre-Relax"},
        copy_files=copy_files,
    )


@job
def mp_metagga_relax_job(
    atoms: Atoms,
    bandgap: float | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Function to relax a structure with Materials Project r2SCAN workflow settings. By default, this uses
    an r2SCAN relax step.

    Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801

    Parameters
    ----------
    atoms
        Atoms object
    bandgap
        Estimate for the bandgap in eV.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VaspSchema
        Dictionary of results.
    """

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
        report_mp_corrections=True,
        additional_fields={"name": "MP Meta-GGA Relax"},
        copy_files=copy_files,
    )


@job
def mp_metagga_static_job(
    atoms: Atoms,
    bandgap: float | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        Files to copy (and decompress) from source to the runtime directory.
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
        report_mp_corrections=True,
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

    1. MP-compatible relax
        - name: "mp_gga_relax_job"
        - job: [quacc.recipes.vasp.mp.mp_gga_relax_job][]

    2. MP-compatible (second) relax
        - name: "mp_gga_relax_job"
        - job: [quacc.recipes.vasp.mp.mp_gga_relax_job][]

    3. MP-compatible static
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

    # Run the second relax
    double_relax_results = mp_gga_relax_job_(
        relax_results["atoms"],
        copy_files={relax_results["dir_name"]: ["CHGCAR*", "WAVECAR*"]},
    )

    # Run the static
    static_results = mp_gga_static_job_(
        double_relax_results["atoms"],
        bandgap=double_relax_results["output"]["bandgap"],
        copy_files={double_relax_results["dir_name"]: ["CHGCAR*", "WAVECAR*"]},
    )

    return {
        "relax1": relax_results,
        "relax2": double_relax_results,
        "static": static_results,
    }


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

    2. MP-compatible relax
        - name: "mp_metagga_relax_job"
        - job: [quacc.recipes.vasp.mp.mp_metagga_relax_job][]

    3. MP-compatible (second) relax
        - name: "mp_metagga_relax_job"
        - job: [quacc.recipes.vasp.mp.mp_metagga_relax_job][]

    4. MP-compatible static
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
    (mp_metagga_prerelax_job_, mp_metagga_relax_job_, mp_metagga_static_job_) = (
        customize_funcs(
            [
                "mp_metagga_prerelax_job",
                "mp_metagga_relax_job",
                "mp_metagga_static_job",
            ],
            [mp_metagga_prerelax_job, mp_metagga_relax_job, mp_metagga_static_job],
            parameters=job_params,
            decorators=job_decorators,
        )
    )

    # Run the prerelax
    prerelax_results = mp_metagga_prerelax_job_(atoms)

    # Run the relax
    relax_results = mp_metagga_relax_job_(
        prerelax_results["atoms"],
        bandgap=prerelax_results["output"]["bandgap"],
        copy_files={prerelax_results["dir_name"]: ["CHGCAR*", "WAVECAR*"]},
    )

    # Run the second relax
    double_relax_results = mp_metagga_relax_job_(
        relax_results["atoms"],
        bandgap=relax_results["output"]["bandgap"],
        copy_files={relax_results["dir_name"]: ["CHGCAR*", "WAVECAR*"]},
    )

    # Run the static
    static_results = mp_metagga_static_job_(
        double_relax_results["atoms"],
        bandgap=double_relax_results["output"]["bandgap"],
        copy_files={double_relax_results["dir_name"]: ["CHGCAR*", "WAVECAR*"]},
    )

    return {
        "prerelax": prerelax_results,
        "relax1": relax_results,
        "relax2": double_relax_results,
        "static": static_results,
    }
