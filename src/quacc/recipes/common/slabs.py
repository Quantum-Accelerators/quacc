"""Common slab workflows."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import subflow
from quacc.atoms.slabs import make_adsorbate_structures, make_slabs_from_bulk

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc import Job


@subflow
def bulk_to_slabs_subflow(
    atoms: Atoms,
    relax_job: Job,
    static_job: Job | None = None,
    make_slabs_kwargs: dict[str, Any] | None = None,
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
    relax_job
        The relaxation function.
    static_job
        The static function.
    make_slabs_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.slabs.make_slabs_from_bulk][]

    Returns
    -------
    list[dict]
        List of schemas.
    """
    make_slabs_kwargs = make_slabs_kwargs or {}

    slabs = make_slabs_from_bulk(atoms, **make_slabs_kwargs)

    results = []
    for slab in slabs:
        result = relax_job(slab)

        if static_job is not None:
            result = static_job(result["atoms"])

        results.append(result)

    return results


@subflow
def slab_to_ads_subflow(
    atoms: Atoms,
    adsorbate: Atoms,
    relax_job: Job,
    static_job: Job | None,
    make_ads_kwargs: dict[str, Any] | None = None,
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
    relax_job
        The slab releaxation job.
    static_job
        The slab static job.
    make_ads_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.slabs.make_adsorbate_structures][]

    Returns
    -------
    list[dict]
        List of schemas.
    """
    make_ads_kwargs = make_ads_kwargs or {}

    slabs = make_adsorbate_structures(atoms, adsorbate, **make_ads_kwargs)

    results = []
    for slab in slabs:
        result = relax_job(slab)

        if static_job is not None:
            result = static_job(result["atoms"])

        results.append(result)

    return results
