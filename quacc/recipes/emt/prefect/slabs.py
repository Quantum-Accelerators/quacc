"""Slab recipes for EMT"""
from __future__ import annotations

from ase import Atoms
from prefect import Flow, Task, flow, task

from quacc.recipes.emt.core import relax_job, static_job
from quacc.util.slabs import make_max_slabs_from_bulk


@flow
def bulk_to_slabs_flow(
    atoms: Atoms,
    slabgen_kwargs: dict | None = None,
    slab_relax_task: Task = task(relax_job),
    slab_static_task: Task | None = task(static_job),
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> Flow:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    slabgen_kwargs
        Additional keyword arguments to pass to make_max_slabs_from_bulk()
    slab_relax_task
        Default Task to use for the relaxation of the slab structures.
    slab_static_task
        Default Task to use for the static calculation of the slab structures.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    Flow
        A Prefect Flow
    """

    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    slabgen_kwargs = slabgen_kwargs or {}

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    def _relax_distributed(slabs):
        return [
            slab_relax_task.submit(slab, **slab_relax_kwargs).result() for slab in slabs
        ]

    def _relax_and_static_distributed(slabs):
        return [
            slab_static_task.submit(
                slab_relax_task.submit(slab, **slab_relax_kwargs).result()["atoms"],
                **slab_static_kwargs,
            ).result()
            for slab in slabs
        ]

    slabs = make_max_slabs_from_bulk(atoms, **slabgen_kwargs)

    if slab_relax_task and slab_static_task:
        return _relax_and_static_distributed(slabs)
    else:
        return _relax_distributed(slabs)
