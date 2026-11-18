"""
Materials Project-compatible recipes using the original MP settings.

This module contains two families of recipes:

- The PBE(+U) recipes (`mp_gga_*`), which use pymatgen's `MPRelaxSet` and
  `MPStaticSet` and the original (unversioned) PBE PAW potentials, which are
  expected in `potpaw_PBE.original` within your `VASP_PP_PATH`.

- The r2SCAN recipes (`mp_prerelax_job`, `mp_metagga_*`), which use pymatgen's
  `MPScanRelaxSet` and `MPScanStaticSet` and the PBE.54 PAW potentials, which are
  expected in `potpaw_PBE.54` within your `VASP_PP_PATH`. These are the settings
  used to compute the r2SCAN data currently in the Materials Project database.
  For the newer MP24 settings, see [quacc.recipes.vasp.mp24][].
"""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires
from pymatgen.io.vasp.sets import MPRelaxSet

from quacc import flow, job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import get_mp_kspacing_handlers, run_and_summarize
from quacc.wflow_tools.customizers import customize_jobs

has_atomate2 = bool(find_spec("atomate2"))

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, TypedDict

    from ase.atoms import Atoms

    from quacc.types import CopyFiles, SourceDirectory, VaspSchema

    class MPGGARelaxFlowSchema(TypedDict):
        """Type hint associated with the MP GGA relaxation flows."""

        relax1: VaspSchema
        relax2: VaspSchema
        static: VaspSchema

    class MPMetaGGARelaxFlowSchema(TypedDict):
        """Type hint associated with the MP meta-GGA relaxation flows."""

        prerelax: VaspSchema
        relax: VaspSchema
        static: VaspSchema


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[atomate2]`")
def mp_gga_relax_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs: Any
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
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results.
    """
    from atomate2.vasp.jobs.mp import MPGGARelaxMaker

    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        MPGGARelaxMaker()
    )
    calc_defaults["incar_copilot_mode"] = "critical"
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        report_mp_corrections=True,
        additional_fields={"name": "MP GGA Relax"},
        copy_files=[{"source": prev_dir, "filenames": ["CHGCAR*", "WAVECAR*"]}]
        if prev_dir
        else None,
    )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[atomate2]`")
def mp_gga_static_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs: Any
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
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
    """
    from atomate2.vasp.jobs.mp import MPGGAStaticMaker

    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        MPGGAStaticMaker()
    )
    calc_defaults["incar_copilot_mode"] = "critical"
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        report_mp_corrections=True,
        additional_fields={"name": "MP GGA Static"},
        copy_files=[{"source": prev_dir, "filenames": ["CHGCAR*", "WAVECAR*"]}]
        if prev_dir
        else None,
    )


@job
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[atomate2]`")
def mp_prerelax_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs: Any
) -> VaspSchema:
    """
    Function to pre-relax a structure with the legacy Materials Project r2SCAN
    workflow settings. By default, this uses a PBEsol pre-relax step.

    These settings match the r2SCAN data currently in the Materials Project
    database and use the PBE.54 PAW potentials.

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
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from atomate2.vasp.jobs.mp import MPPreRelaxMaker

    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        MPPreRelaxMaker()
    )
    calc_defaults["incar_copilot_mode"] = "critical"
    calc_defaults["custodian_handlers"] = get_mp_kspacing_handlers()
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP Legacy PBEsol Pre-Relax"},
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
    Function to relax a structure with the legacy Materials Project r2SCAN
    workflow settings. By default, this uses an r2SCAN relax step.

    These settings match the r2SCAN data currently in the Materials Project
    database and use the PBE.54 PAW potentials.

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
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from atomate2.vasp.jobs.mp import MPMetaGGARelaxMaker

    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        MPMetaGGARelaxMaker()
    )
    calc_defaults["incar_copilot_mode"] = "critical"
    calc_defaults["custodian_handlers"] = get_mp_kspacing_handlers()
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP Legacy r2SCAN Relax"},
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
    Function to run a static calculation on a structure with the legacy Materials
    Project r2SCAN workflow settings. By default, this uses an r2SCAN static step.

    These settings match the r2SCAN data currently in the Materials Project
    database and use the PBE.54 PAW potentials.

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
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][]. All of the ASE
        Vasp calculator keyword arguments are supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
        See the type-hint for the data structure.
    """
    from atomate2.vasp.jobs.mp import MPMetaGGAStaticMaker

    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_maker(
        MPMetaGGAStaticMaker()
    )
    calc_defaults["incar_copilot_mode"] = "critical"
    calc_defaults["custodian_handlers"] = get_mp_kspacing_handlers()
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP Legacy r2SCAN Static"},
        copy_files=[{"source": prev_dir, "filenames": ["CHGCAR*", "WAVECAR*"]}]
        if prev_dir
        else None,
    )


@job
def mp_relax_set_job(
    atoms: Atoms,
    copy_files: CopyFiles | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
) -> VaspSchema:
    """
    Run a static calculation with the original pymatgen `MPRelaxSet` VASP settings,
    compatible with the MPtrj, WBM, sAlex, and MatterSim datasets. As in the
    Materials Project relaxation workflows, Gamma-centered k-point grids are used.

    Parameters
    ----------
    atoms
        Atoms object.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Vasp calculator. Set a value to `None` to remove a
        pre-existing key entirely. All ASE Vasp calculator keyword arguments are
        supported.

    Returns
    -------
    VaspSchema
        Dictionary of results from [quacc.schemas.vasp.VaspSummarize.run][].
    """
    # Materials Project relaxations (atomate and atomate2) use Gamma-centered grids
    calc_defaults = MPtoASEConverter(atoms=atoms).convert_input_set(
        MPRelaxSet(force_gamma=True)
    )
    calc_defaults |= {"nsw": 0, "incar_copilot_mode": "critical"}

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MP RelaxSet Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@flow
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[atomate2]`")
def mp_gga_relax_flow(
    atoms: Atoms,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> MPGGARelaxFlowSchema:
    """
    Materials Project GGA workflow consisting of:

    1. MP-compatible relax
        - name: "mp_gga_relax_job"
        - job: [quacc.recipes.vasp.mp_legacy.mp_gga_relax_job][]

    2. MP-compatible (second) relax
        - name: "mp_gga_relax_job"
        - job: [quacc.recipes.vasp.mp_legacy.mp_gga_relax_job][]

    3. MP-compatible static
        - name: "mp_gga_static_job"
        - job: [quacc.recipes.vasp.mp_legacy.mp_gga_static_job][]

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
    (mp_gga_relax_job_, mp_gga_static_job_) = customize_jobs(
        {"mp_gga_relax_job": mp_gga_relax_job, "mp_gga_static_job": mp_gga_static_job},
        param_swaps=job_params,
        decorators=job_decorators,
    ).values()

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
@requires(has_atomate2, "atomate2 is not installed. Run `pip install quacc[atomate2]`")
def mp_metagga_relax_flow(
    atoms: Atoms,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> MPMetaGGARelaxFlowSchema:
    """
    Legacy Materials Project r2SCAN workflow, matching atomate2's
    `MPMetaGGADoubleRelaxStaticMaker`. These settings match the r2SCAN data
    currently in the Materials Project database. It consists of:

    1. MP-compatible PBEsol pre-relax
        - name: "mp_prerelax_job"
        - job: [quacc.recipes.vasp.mp_legacy.mp_prerelax_job][]

    2. MP-compatible r2SCAN relax
        - name: "mp_metagga_relax_job"
        - job: [quacc.recipes.vasp.mp_legacy.mp_metagga_relax_job][]

    3. MP-compatible r2SCAN static
        - name: "mp_metagga_static_job"
        - job: [quacc.recipes.vasp.mp_legacy.mp_metagga_static_job][]

    Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801

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
