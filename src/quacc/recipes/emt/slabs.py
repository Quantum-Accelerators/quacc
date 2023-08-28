"""Slab recipes for EMT"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow, job, subflow
from quacc.recipes.emt.core import relax_job as _relax_job
from quacc.recipes.emt.core import static_job as _static_job
from quacc.utils.slabs import make_max_slabs_from_bulk
from quacc.utils.wflows import fetch_atoms

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import OptSchema, RunSchema
    from quacc.utils.wflows import Job


@flow
def bulk_to_slabs_flow(
    atoms: Atoms | dict,
    make_slabs_kwargs: dict | None = None,
    slab_relax: Job = _relax_job,
    slab_static: Job | None = _static_job,
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> list[RunSchema | OptSchema]:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    make_slabs_kwargs
        Additional keyword arguments to pass to `quacc.utils.slabs.make_max_slabs_from_bulk`
    slab_relax
        Default Job to use for the relaxation of the slab structures.
    slab_static
        Default Job to use for the static calculation of the slab structures.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    list[RunSchema | OptSchema]
        RunSchema or OptSchema for each slab.
    """
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    make_slabs_kwargs = make_slabs_kwargs or {}

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    @job
    def _make_slabs(atoms):
        atoms = fetch_atoms(atoms)
        return make_max_slabs_from_bulk(atoms, **make_slabs_kwargs)

    @subflow
    def _relax_distributed(slabs):
        return [slab_relax(slab, **slab_relax_kwargs) for slab in slabs]

    @subflow
    def _relax_and_static_distributed(slabs):
        return [
            slab_static(
                slab_relax(slab, **slab_relax_kwargs),
                **slab_static_kwargs,
            )
            for slab in slabs
        ]

    slabs = _make_slabs(atoms)

    if slab_static is None:
        return _relax_distributed(slabs)

    return _relax_and_static_distributed(slabs)
