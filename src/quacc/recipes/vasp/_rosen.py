"""
Internal Rosen recipes in development.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.io.vasp.sets import MatPESStaticSet

from quacc import flow, job
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.recipes.vasp._base import run_and_summarize
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, TypedDict

    from ase.atoms import Atoms

    from quacc.types import SourceDirectory, VaspSchema

    class MLIPFlowSchema(TypedDict):
        gga: VaspSchema
        metagga: VaspSchema | None


@job
def gga_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    GGA single point calculation.

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
    dict_set = MatPESStaticSet(
        auto_kspacing=False,
        auto_ispin=True,
        user_incar_settings={
            "ALGO": "All",
            "EFERMI": "Midgap",
            "ENAUG": None,
            "GGA_COMPAT": False,
            "KSPACING": 0.4,
            "LELF": True,
            "LWAVE": True,
        },
    )
    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_dict_set(
        dict_set
    )

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "GGA SP"},
        copy_files={prev_dir: ["WAVECAR*"]},
    )


@job
def metagga_job(
    atoms: Atoms, prev_dir: SourceDirectory | None = None, **calc_kwargs
) -> VaspSchema:
    """
    MetaGGA single point calculation.


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
    dict_set = MatPESStaticSet(
        xc_functional="R2SCAN",
        auto_kspacing=False,
        auto_ispin=True,
        user_incar_settings={
            "ALGO": "All",
            "EFERMI": "Midgap",
            "ENAUG": None,
            "KSPACING": 0.4,
            "GGA_COMPAT": False,
            "LELF": True,
        },
    )
    calc_defaults = MPtoASEConverter(atoms=atoms, prev_dir=prev_dir).convert_dict_set(
        dict_set
    )
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "MetaGGA SP"},
        copy_files={prev_dir: ["WAVECAR*"]},
    )


@flow
def mlip_flow(
    atoms: Atoms,
    run_metagga: bool = True,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> MLIPFlowSchema:
    """
    Materials Project r2SCAN workflow consisting of:

    1. GGA SP
        - name: "gga_job"
        - job: [quacc.recipes.vasp.mp.mp_prerelax_job][]

    2. MP-compatible relax
        - name: "metagga_job"
        - job: [quacc.recipes.vasp.mp.mp_metagga_relax_job][]

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    run_metagga
        Whether to run the MetaGGA job.
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
    (gga_job_, metagga_job_) = customize_funcs(
        ["gga_job", "metagga_job"],
        [gga_job, metagga_job],
        param_swaps=job_params,
        decorators=job_decorators,
    )

    gga_results = gga_job_(atoms)
    if run_metagga:
        metagga_results = metagga_job_(atoms, prev_dir=gga_results["dir_name"])
    else:
        metagga_results = None

    return {"gga": gga_results, "metagga": metagga_results}
