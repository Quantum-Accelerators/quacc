"""Recipes for slabs."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow, job
from quacc.recipes.common.slabs import bulk_to_slabs_subflow, slab_to_ads_subflow
from quacc.recipes.vasp._base import base_fn
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.schemas._aliases.vasp import VaspSchema
    from quacc.utils.files import Filenames, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    preset: str | None = "SlabSet",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Function to carry out a single-point calculation on a slab.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
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
        "auto_dipole": True,
        "ismear": -5,
        "laechg": True,
        "lcharg": True,
        "lreal": False,
        "lvhar": True,
        "lwave": True,
        "nedos": 3001,
        "nsw": 0,
    }
    return base_fn(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "VASP Slab Static"},
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    preset: str | None = "SlabSet",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> VaspSchema:
    """
    Function to relax a slab.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
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
        "auto_dipole": True,
        "ediffg": -0.02,
        "isif": 2,
        "ibrion": 2,
        "isym": 0,
        "lcharg": False,
        "lwave": False,
        "nsw": 200,
        "symprec": 1e-8,
    }
    return base_fn(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "VASP Slab Relax"},
        copy_files=copy_files,
    )


@flow
def bulk_to_slabs_flow(
    atoms: Atoms,
    make_slabs_kwargs: dict[str, Any] | None = None,
    run_static: bool = True,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> list[VaspSchema]:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations
        - name: "relax_job"
        - job: [quacc.recipes.vasp.slabs.relax_job][]

    3. Optional slab statics
        - name: "static_job"
        - job: [quacc.recipes.vasp.slabs.static_job][]

    Parameters
    ----------
    atoms
        Atoms object
    make_slabs_kwargs
        Additional keyword arguments to pass to [quacc.atoms.slabs.make_slabs_from_bulk][]
    run_static
        Whether to run static calculations.
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictinoary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    list[VaspSchema]
        List of dictionary results from [quacc.schemas.vasp.vasp_summarize_run][].
        See the type-hint for the data structure.
    """
    relax_job_, static_job_ = customize_funcs(
        ["relax_job", "static_job"],
        [relax_job, static_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    return bulk_to_slabs_subflow(
        atoms,
        relax_job_,
        static_job=static_job_ if run_static else None,
        make_slabs_kwargs=make_slabs_kwargs,
    )


@flow
def slab_to_ads_flow(
    slab: Atoms,
    adsorbate: Atoms,
    run_static: bool = True,
    make_ads_kwargs: dict[str, Any] | None = None,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> list[VaspSchema]:
    """
    Workflow consisting of:

    1. Slab-adsorbate generation

    2. Slab-adsorbate relaxations
        - name: "relax_job"
        - job: [quacc.recipes.vasp.slabs.relax_job][]

    3. Optional slab-adsorbate statics
        - name: "static_job"
        - job: [quacc.recipes.vasp.slabs.static_job][]

    Parameters
    ----------
    slab
        Atoms object for the slab structure.
    adsorbate
        Atoms object for the adsorbate.
    run_static
        Whether to run static calculations.
    make_ads_kwargs
        Additional keyword arguments to pass to [quacc.atoms.slabs.make_adsorbate_structures][]
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictinoary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    list[VaspSchema]
        List of dictionaries of results from [quacc.schemas.vasp.vasp_summarize_run][].
        See the type-hint for the data structure.
    """
    relax_job_, static_job_ = customize_funcs(
        ["relax_job", "static_job"],
        [relax_job, static_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    return slab_to_ads_subflow(
        slab,
        adsorbate,
        relax_job_,
        static_job=static_job_ if run_static else None,
        make_ads_kwargs=make_ads_kwargs,
    )
