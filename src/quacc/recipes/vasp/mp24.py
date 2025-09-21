"""
Materials Project-compatible recipes using the MP24 sets.

!!! Important

    Make sure that you use the Materials Project-compatible pseudpotential
    versions (i.e. v.64)
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
    from atomate2.vasp.jobs.mp import MP24PreRelaxMaker, MP24RelaxMaker, MP24StaticMaker

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, TypedDict

    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema

    class MPMetaGGARelaxFlowSchema(TypedDict):
        """Type hint associated with the MP meta-GGA relaxation flows."""

        prerelax: VaspSchema
        relax1: VaspSchema
        relax2: VaspSchema
        static: VaspSchema


_MP_SETTINGS = {"VASP_INCAR_COPILOT": "off", "VASP_USE_CUSTODIAN": True}


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mp_prerelax_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Function to pre-relax a structure with Materials Project r2SCAN workflow settings. By default, this
    uses a PBEsol pre-relax step.

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
        MP24PreRelaxMaker()
    )
    with change_settings(_MP_SETTINGS):
        return run_and_summarize(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            report_mp_corrections=True,
            additional_fields={"name": "MP PBESol Pre-Relax"},
            copy_files={prev_dir: ["WAVECAR*"]} if prev_dir else None,
        )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[mp]`")
def mp_metagga_relax_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    Function to relax a structure with Materials Project r2SCAN workflow settings. By default, this uses
    an r2SCAN relax step.

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
        MP24RelaxMaker()
    )
    with change_settings(_MP_SETTINGS):
        return run_and_summarize(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            report_mp_corrections=True,
            additional_fields={"name": "MP r2SCAN Relax"},
            copy_files={prev_dir: ["WAVECAR*"]} if prev_dir else None,
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
        MP24StaticMaker()
    )
    with change_settings(_MP_SETTINGS):
        return run_and_summarize(
            atoms,
            calc_defaults=calc_defaults,
            calc_swaps=calc_kwargs,
            report_mp_corrections=True,
            additional_fields={"name": "MP r2SCAN Static"},
            copy_files={prev_dir: ["WAVECAR*"]} if prev_dir else None,
        )


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
        - name: "mpa_prerelax_job"
        - job: [quacc.recipes.vasp.mp24.mp_prerelax_job][]

    2. MP-compatible relax
        - name: "mp_metagga_relax_job"
        - job: [quacc.recipes.vasp.mp24.mp_metagga_relax_job][]

    3. MP-compatible static
        - name: "mp_metagga_static_job"
        - job: [quacc.recipes.vasp.mp24.mp_metagga_static_job][]

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
    (mp_prerelax_job_, mp_metagga_relax_job_, mp_metagga_static_job_) = customize_funcs(
        ["mp_prerelax_job", "mp_metagga_relax_job", "mp_metagga_static_job"],
        [mp_prerelax_job, mp_metagga_relax_job, mp_metagga_static_job],
        param_swaps=job_params,
        decorators=job_decorators,
    )

    # Run the prerelax
    prerelax_results = mp_prerelax_job_(atoms)

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
