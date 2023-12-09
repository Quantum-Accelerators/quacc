"""Recipes for slabs."""
from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

from quacc import flow, job
from quacc.atoms.slabs import make_adsorbate_structures, make_slabs_from_bulk
from quacc.recipes.common.slabs import bulk_to_slabs_subflow, slab_to_ads_subflow
from quacc.recipes.vasp._base import base_fn

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc import Job
    from quacc.schemas._aliases.vasp import VaspSchema


@job
def static_job(
    atoms: Atoms,
    preset: str | None = "SlabSet",
    copy_files: list[str] | None = None,
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
        "auto_dipole": True,
        "ismear": -5,
        "laechg": True,
        "lcharg": True,
        "lreal": False,
        "lvhar": True,
        "lwave": True,
        "nedos": 5001,
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
    copy_files: list[str] | None = None,
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
    slab_relax_job: Job = relax_job,
    slab_static_job: Job | None = static_job,
) -> list[VaspSchema]:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    atoms
        Atoms object
    make_slabs_kwargs
        Additional keyword arguments to pass to [quacc.atoms.slabs.make_slabs_from_bulk][]
    slab_relax_job
        Relaxation job, which defaults to [quacc.recipes.vasp.slabs.slab_relax_job][].
    slab_static_kwargs
        Static job, which defaults to [quacc.recipes.vasp.slabs.slab_static_job][].

    Returns
    -------
    list[VaspSchema]
        List of dictionary results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    make_slabs_kwargs = make_slabs_kwargs or {}
    return bulk_to_slabs_subflow(
        atoms,
        slab_relax_job,
        static_job=slab_static_job,
        make_slabs_fn=partial(make_slabs_from_bulk, **make_slabs_kwargs),
    )


@flow
def slab_to_ads_flow(
    slab: Atoms,
    adsorbate: Atoms,
    make_ads_kwargs: dict[str, Any] | None = None,
    slab_relax_job: Job = relax_job,
    slab_static_job: Job | None = static_job,
) -> list[VaspSchema]:
    """
    Workflow consisting of:

    1. Slab-adsorbate generation

    2. Slab-adsorbate relaxations

    3. Slab-adsorbate statics (optional)

    Parameters
    ----------
    slab
        Atoms object for the slab structure.
    adsorbate
        Atoms object for the adsorbate.
    make_ads_kwargs
        Additional keyword arguments to pass to [quacc.atoms.slabs.make_adsorbate_structures][]
    slab_relax_kwargs
        Relaxation job, which defaults to [quacc.recipes.vasp.slabs.slab_relax_job][].
    slab_static_kwargs
        Static job, which defaults to [quacc.recipes.vasp.slabs.slab_static_job][].

    Returns
    -------
    list[VaspSchema]
        List of dictionaries of results from [quacc.schemas.vasp.vasp_summarize_run][]
    """

    make_ads_kwargs = make_ads_kwargs or {}
    return slab_to_ads_subflow(
        slab,
        adsorbate,
        slab_relax_job,
        static_job=slab_static_job,
        make_ads_fn=partial(make_adsorbate_structures, **make_ads_kwargs),
    )
