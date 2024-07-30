"""
Materials Project-compatible recipes.

!!! Important

    Make sure that you use the Materials Project-compatible pseudpotential
    versions. The GGA workflows use the old (no version) PAW PBE potentials.
    The meta-GGA workflows currently use the v.54 PAW PBE potentials.
"""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import change_settings, flow, job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import run_and_summarize
from quacc.wflow_tools.customizers import customize_funcs

has_atomate2 = bool(find_spec("atomate2"))
if has_atomate2:
    from atomate2.vasp.jobs.mp import (
        MPGGARelaxMaker,
        MPGGAStaticMaker,
        MPMetaGGARelaxMaker,
        MPMetaGGAStaticMaker,
        MPPreRelaxMaker,
    )

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.types import (
        MPGGARelaxFlowSchema,
        MPMetaGGARelaxFlowSchema,
        SourceDirectory,
        VaspSchema,
    )

_MP_SETTINGS = {"VASP_INCAR_COPILOT": "off", "VASP_USE_CUSTODIAN": True}


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mp_gga_relax_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Function to relax a structure with the original Materials Project GGA(+U) settings.

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        A previous directory for a prior step in the workflow.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VaspSchema
        Dictionary of results.
    """
    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_vasp_maker(
        MPGGARelaxMaker
    )
    with change_settings(_MP_SETTINGS):
        return run_and_summarize(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            report_mp_corrections=True,
            additional_fields={"name": "MP GGA Relax"},
            copy_files={prev_dir: ["CHGCAR*", "WAVECAR*"]},
        )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mp_gga_static_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Function to run a static calculation on a structure with the original Materials Project GGA(+U) settings.

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        A previous directory for a prior step in the workflow.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
    """
    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_vasp_maker(
        MPGGAStaticMaker
    )
    with change_settings(_MP_SETTINGS):
        return run_and_summarize(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            report_mp_corrections=True,
            additional_fields={"name": "MP GGA Static"},
            copy_files={prev_dir: ["CHGCAR*", "WAVECAR*"]},
        )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mp_metagga_prerelax_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Function to pre-relax a structure with Materials Project r2SCAN workflow settings. By default, this
    uses a PBEsol pre-relax step.

    Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        A previous directory for a prior step in the workflow.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_vasp_maker(
        MPPreRelaxMaker
    )
    with change_settings(_MP_SETTINGS):
        return run_and_summarize(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            report_mp_corrections=True,
            additional_fields={"name": "MP Meta-GGA Pre-Relax"},
            copy_files={prev_dir: ["CHGCAR*", "WAVECAR*"]},
        )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mp_metagga_relax_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Function to relax a structure with Materials Project r2SCAN workflow settings. By default, this uses
    an r2SCAN relax step.

    Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        A previous directory for a prior step in the workflow.
    **calc_kwargs
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].

    Returns
    -------
    VaspSchema
        Dictionary of results.
    """
    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_vasp_maker(
        MPMetaGGARelaxMaker
    )
    with change_settings(_MP_SETTINGS):
        return run_and_summarize(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            report_mp_corrections=True,
            additional_fields={"name": "MP Meta-GGA Relax"},
            copy_files={prev_dir: ["CHGCAR*", "WAVECAR*"]},
        )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mp_metagga_static_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Function to run a static calculation on a structure with r2SCAN workflow Materials Project settings.
    By default, this uses an r2SCAN static step.

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        A previous directory for a prior step in the workflow.
    **calc_kwargs
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to `ase.calculators.vasp.vasp.Vasp`.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_vasp_maker(
        MPMetaGGAStaticMaker
    )
    with change_settings(_MP_SETTINGS):
        return run_and_summarize(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            report_mp_corrections=True,
            additional_fields={"name": "MP Meta-GGA Static"},
            copy_files={prev_dir: ["CHGCAR*", "WAVECAR*"]},
        )


@flow
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
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
        param_swaps=job_params,
        decorators=job_decorators,
    )

    # Run the relax
    relax_results = mp_gga_relax_job_(atoms)

    # Run the second relax
    double_relax_results = mp_gga_relax_job_(
        relax_results["atoms"], prev_dir=relax_results["dir_name"]
    )

    # Run the static
    static_results = mp_gga_static_job_(
        double_relax_results["atoms"], prev_dir=double_relax_results["dir_name"]
    )

    return {
        "relax1": relax_results,
        "relax2": double_relax_results,
        "static": static_results,
    }


@flow
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
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
            param_swaps=job_params,
            decorators=job_decorators,
        )
    )

    # Run the prerelax
    prerelax_results = mp_metagga_prerelax_job_(atoms)

    # Run the relax
    relax_results = mp_metagga_relax_job_(
        prerelax_results["atoms"], prev_dir=prerelax_results["dir_name"]
    )

    # Run the second relax
    double_relax_results = mp_metagga_relax_job_(
        relax_results["atoms"], prev_dir=relax_results["dir_name"]
    )

    # Run the static
    static_results = mp_metagga_static_job_(
        double_relax_results["atoms"], prev_dir=double_relax_results["dir_name"]
    )

    return {
        "prerelax": prerelax_results,
        "relax1": relax_results,
        "relax2": double_relax_results,
        "static": static_results,
    }
