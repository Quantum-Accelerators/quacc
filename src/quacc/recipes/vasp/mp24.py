"""
Materials Project-compatible recipes using the MP24 sets.

These recipes use the r2SCAN settings that the Materials Project adopted after its
2023-2024 benchmarking effort (`MP24RelaxSet` and `MP24StaticSet` in pymatgen),
which use the PBE.64 PAW potentials. To generate calculations that are compatible
with the r2SCAN data currently in the Materials Project database, which was
computed with the PBE.54 PAW potentials, use the r2SCAN recipes in
[quacc.recipes.vasp.mp_legacy][] instead.

!!! Important

    Your `VASP_PP_PATH` must contain the `potpaw_PBE.64` pseudopotential library.
"""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import flow, job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import get_mp_kspacing_handlers, run_and_summarize
from quacc.wflow_tools.customizers import customize_jobs

has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, TypedDict

    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema

    class MPMetaGGARelaxFlowSchema(TypedDict):
        """Type hint associated with the MP meta-GGA relaxation flows."""

        prerelax: VaspSchema
        relax: VaspSchema
        static: VaspSchema


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[atomate2]`")
def mp_prerelax_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs: Any
) -> VaspSchema:
    """
    Function to pre-relax a structure with the Materials Project MP24 settings. By
    default, this uses a PBEsol pre-relax step.

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        A previous directory for a prior step in the workflow.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from atomate2.vasp.jobs.mp import MP24PreRelaxMaker

    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        MP24PreRelaxMaker()
    )
    calc_defaults["incar_copilot_mode"] = "critical"
    calc_defaults["custodian_handlers"] = get_mp_kspacing_handlers()
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP PBESol Pre-Relax"},
        copy_files=[{"source": prev_dir, "filenames": ["CHGCAR*", "WAVECAR*"]}]
        if prev_dir
        else None,
    )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[atomate2]`")
def mp_metagga_relax_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs: Any
) -> VaspSchema:
    """
    Function to relax a structure with the Materials Project MP24 settings. By
    default, this uses an r2SCAN relax step.

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        A previous directory for a prior step in the workflow.
    **calc_kwargs
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results.
    """
    from atomate2.vasp.jobs.mp import MP24RelaxMaker

    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        MP24RelaxMaker()
    )
    calc_defaults["incar_copilot_mode"] = "critical"
    calc_defaults["custodian_handlers"] = get_mp_kspacing_handlers()
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP r2SCAN Relax"},
        copy_files=[{"source": prev_dir, "filenames": ["CHGCAR*", "WAVECAR*"]}]
        if prev_dir
        else None,
    )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[atomate2]`")
def mp_metagga_static_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs: Any
) -> VaspSchema:
    """
    Function to run a static calculation on a structure with the Materials Project
    MP24 settings. By default, this uses an r2SCAN static step.

    Parameters
    ----------
    atoms
        Atoms object
    prev_dir
        A previous directory for a prior step in the workflow.
    **calc_kwargs
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from atomate2.vasp.jobs.mp import MP24StaticMaker

    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        MP24StaticMaker()
    )
    calc_defaults["incar_copilot_mode"] = "critical"
    calc_defaults["custodian_handlers"] = get_mp_kspacing_handlers()
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP r2SCAN Static"},
        copy_files=[{"source": prev_dir, "filenames": ["CHGCAR*", "WAVECAR*"]}]
        if prev_dir
        else None,
    )


@flow
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[atomate2]`")
def mp_metagga_relax_flow(
    atoms: Atoms,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> MPMetaGGARelaxFlowSchema:
    """
    Materials Project MP24 r2SCAN workflow, matching atomate2's
    `MP24DoubleRelaxStaticMaker`. It consists of:

    1. MP24 PBEsol pre-relax
        - name: "mp_prerelax_job"
        - job: [quacc.recipes.vasp.mp24.mp_prerelax_job][]

    2. MP24 r2SCAN relax
        - name: "mp_metagga_relax_job"
        - job: [quacc.recipes.vasp.mp24.mp_metagga_relax_job][]

    3. MP24 r2SCAN static
        - name: "mp_metagga_static_job"
        - job: [quacc.recipes.vasp.mp24.mp_metagga_static_job][]

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    MPMetaGGARelaxFlowSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    (mp_prerelax_job_, mp_metagga_relax_job_, mp_metagga_static_job_) = customize_jobs(
        {
            "mp_prerelax_job": mp_prerelax_job,
            "mp_metagga_relax_job": mp_metagga_relax_job,
            "mp_metagga_static_job": mp_metagga_static_job,
        },
        param_swaps=job_params,
        decorators=job_decorators,
    ).values()

    # Run the PBEsol pre-relax
    prerelax_results = mp_prerelax_job_(atoms)

    # Run the r2SCAN relax
    relax_results = mp_metagga_relax_job_(
        prerelax_results["atoms"], prev_dir=prerelax_results["dir_name"]
    )

    # Run the r2SCAN static
    static_results = mp_metagga_static_job_(
        relax_results["atoms"], prev_dir=relax_results["dir_name"]
    )

    return {
        "prerelax": prerelax_results,
        "relax": relax_results,
        "static": static_results,
    }
