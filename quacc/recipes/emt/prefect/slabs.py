"""Slab recipes for EMT"""
from __future__ import annotations

from ase import Atoms
from prefect import flow

from quacc.util.slabs import make_max_slabs_from_bulk


# TODO: Add type hints. Relies on #9266 and #10116 in Prefect
# TODO: Make `slab_relax_task` and `slab_static_task` kwargs
# like in other workflow engines. Relies on #10135 in Prefect
@flow
def bulk_to_slabs_flow(
    atoms,
    slab_relax_task,
    slab_static_task,
    slabgen_kwargs: dict | None = None,
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> list[dict]:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value.
    slab_relax_task
        Default Task to use for the relaxation of the slab structures.
    slab_static_task
        Default Task to use for the static calculation of the slab structures.
    slabgen_kwargs
        Additional keyword arguments to pass to make_max_slabs_from_bulk()
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    list[dict]
        List of dictionary of results from quacc.schemas.ase.summarize_run
        or quacc.schemas.ase.summarize_opt_run
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    slabgen_kwargs = slabgen_kwargs or {}

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    def _relax_distributed(slabs):
        return [slab_relax_task(slab, **slab_relax_kwargs) for slab in slabs]

    def _relax_and_static_distributed(slabs):
        return [
            slab_static_task.submit(
                slab_relax_task.submit(slab, **slab_relax_kwargs).result(),
                **slab_static_kwargs,
            ).result()
            for slab in slabs
        ]

    slabs = make_max_slabs_from_bulk(atoms, **slabgen_kwargs)

    if slab_static_task is None:
        return _relax_distributed(slabs)

    return _relax_and_static_distributed(slabs)
