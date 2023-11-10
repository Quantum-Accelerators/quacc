"""Common slab workflows"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import subflow
from quacc.atoms.slabs import make_slabs_from_bulk

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema


@subflow
def common_bulk_to_slabs_flow(
    atoms: Atoms,
    relax_job: Callable,
    static_job: Callable | None,
    make_slabs_kwargs: dict[str, Any] | None = None,
    relax_job_kwargs: dict[str, Any] | None = None,
    static_job_kwargs: dict[str, Any] | None = None,
) -> list[RunSchema | OptSchema]:
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
    relax_job_kwargs
        Additional keyword arguments to pass to [quacc.recipes.emt.core.relax_job][].
    static_job_kwargs
        Additional keyword arguments to pass to [quacc.recipes.emt.core.static_job][].

    Returns
    -------
    list[RunSchema | OptSchema]
        [RunSchema][quacc.schemas.ase.summarize_run] or
        [OptSchema][quacc.schemas.ase.summarize_opt_run] for each slab.
    """
    relax_job_kwargs = relax_job_kwargs or {}
    static_job_kwargs = static_job_kwargs or {}
    make_slabs_kwargs = make_slabs_kwargs or {}

    slabs = make_slabs_from_bulk(atoms, **make_slabs_kwargs)

    results = []
    for slab in slabs:
        result = relax_job(slab, **relax_job_kwargs)

        if static_job:
            result = static_job(result["atoms"], **static_job_kwargs)

        results.append(result)

    return results
