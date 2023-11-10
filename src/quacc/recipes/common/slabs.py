"""Common slab workflows"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import subflow
from quacc.atoms.slabs import make_adsorbate_structures, make_slabs_from_bulk

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase import Atoms


@subflow
def bulk_to_slabs_subflow(
    atoms: Atoms,
    slab_relax_job: Callable,
    slab_static_job: Callable | None,
    make_slabs_kwargs: dict[str, Any] | None = None,
    slab_relax_kwargs: dict[str, Any] | None = None,
    slab_static_kwargs: dict[str, Any] | None = None,
) -> list[dict]:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    atoms
        Atoms object
    slab_relax_job
        The relaxation function.
    slab_static_job
        The static function.
    make_slabs_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.slabs.make_slabs_from_bulk][]
    slab_relax_kwargs
        Additional keyword arguments to pass to [quacc.recipes.emt.core.relax_job][].
    slab_static_kwargs
        Additional keyword arguments to pass to [quacc.recipes.emt.core.static_job][].

    Returns
    -------
    list[dict]
        List of schemas.
    """
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    make_slabs_kwargs = make_slabs_kwargs or {}

    slabs = make_slabs_from_bulk(atoms, **make_slabs_kwargs)

    results = []
    for slab in slabs:
        result = slab_relax_job(slab, **slab_relax_kwargs)

        if slab_static_job:
            result = slab_static_job(result["atoms"], **slab_static_kwargs)

        results.append(result)

    return results


@subflow
def slab_to_ads_subflow(
    atoms: Atoms,
    adsorbate: Atoms,
    slab_relax_job: Callable,
    slab_static_job: Callable | None,
    make_ads_kwargs: dict[str, Any] | None = None,
    slab_relax_kwargs: dict[str, Any] | None = None,
    slab_static_kwargs: dict[str, Any] | None = None,
) -> list[dict]:
    """
    Workflow consisting of:

    1. Slab-adsorbate generation

    2. Slab-adsorbate relaxations

    3. Slab-adsorbate statics (optional)

    Parameters
    ----------
    atoms
        Atoms object for the slab structure.
    adsorbate
        Atoms object for the adsorbate.
    slab_relax_job
        The slab releaxation job.
    slab_static_job
        The slab static job.
    make_ads_kwargs
        Additional keyword arguments to pass to [quacc.atoms.slabs.make_adsorbate_structures][]
    slab_relax_kwargs
        Additional keyword arguments to pass to [quacc.recipes.vasp.slabs.relax_job][].
    slab_static_kwargs
        Additional keyword arguments to pass to [quacc.recipes.vasp.slabs.static_job][].

    Returns
    -------
    list[dict]
        List of schemas.
    """

    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    make_ads_kwargs = make_ads_kwargs or {}

    slabs = make_adsorbate_structures(atoms, adsorbate, **make_ads_kwargs)

    results = []
    for slab in slabs:
        result = slab_relax_job(slab, **slab_relax_kwargs)

        if slab_static_job:
            result = slab_static_job(result["atoms"], **slab_static_kwargs)

        results.append(result)

    return results
