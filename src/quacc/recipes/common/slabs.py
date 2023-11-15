"""Common slab workflows"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import subflow
from quacc.atoms.slabs import make_adsorbate_structures, make_slabs_from_bulk

if TYPE_CHECKING:
    from typing import Callable

    from ase import Atoms


@subflow
def bulk_to_slabs_subflow(
    atoms: Atoms,
    relax_fn: Callable,
    static_fn: Callable | None = None,
    make_slabs_fn: Callable = make_slabs_from_bulk,
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
    relax_fn
        The relaxation function.
    static_fn
        The static function.
    make_slabs_fn
        The function for generating slabs.

    Returns
    -------
    list[dict]
        List of schemas.
    """

    slabs = make_slabs_fn(atoms)

    results = []
    for slab in slabs:
        result = relax_fn(slab)

        if static_fn:
            result = static_fn(result["atoms"])

        results.append(result)

    return results


@subflow
def slab_to_ads_subflow(
    atoms: Atoms,
    adsorbate: Atoms,
    relax_fn: Callable,
    static_fn: Callable | None,
    make_ads_fn: Callable = make_adsorbate_structures,
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
    relax_fn
        The slab releaxation job.
    static_fn
        The slab static job.
    make_ads_fn
        The function to generate slab-adsorbate structures.

    Returns
    -------
    list[dict]
        List of schemas.
    """

    slabs = make_ads_fn(atoms, adsorbate)

    results = []
    for slab in slabs:
        result = relax_fn(slab)

        if static_fn:
            result = static_fn(result["atoms"])

        results.append(result)

    return results
