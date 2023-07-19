"""Slab recipes for EMT"""
from __future__ import annotations

from ase import Atoms
from prefect import flow, task
from prefect.futures import PrefectFuture

from quacc.recipes.emt.slabs import relax_job, static_job
from quacc.schemas.ase import OptSchema, RunSchema
from quacc.util.slabs import make_max_slabs_from_bulk


@flow
def bulk_to_slabs_flow(
    atoms,
    run_slab_static: bool = True,
    make_slabs_kwargs: dict | None = None,
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> list[PrefectFuture]:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value.
    run_slab_static
        Whether to run the slab static calculations.
    make_slabs_kwargs
        Additional keyword arguments to pass to make_max_slabs_from_bulk()
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    list[PrefectFuture]
        List of PrefectFuture objects, each of which resolves to a dictionary of results
        from quacc.schemas.ase.summarize_run or quacc.schemas.ase.summarize_opt_run
    """
    make_slabs_kwargs = make_slabs_kwargs or {}
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}

    slab_relax_task = task(relax_job)
    slab_static_task = task(static_job)

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    # Generate all the slab
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    slabs = make_max_slabs_from_bulk(atoms, **make_slabs_kwargs)

    futures = []
    for slab in slabs:
        slab_relax_future = slab_relax_task.submit(slab, **slab_relax_kwargs)
        if run_slab_static:
            slab_static_future = slab_static_task.submit(
                slab_relax_future, **slab_static_kwargs
            )
            futures.append(slab_static_future)
        else:
            futures.append(slab_relax_future)

    return futures
