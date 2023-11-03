"""Slab recipes for EMT."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow, job, subflow
from quacc.atoms.slabs import make_slabs_from_bulk
from quacc.recipes.emt.core import relax_job, static_job

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas.ase import OptSchema, RunSchema


@flow
def bulk_to_slabs_flow(
    atoms: Atoms,
    make_slabs_kwargs: dict[str, Any] | None = None,
    run_static: bool = True,
    slab_relax_kwargs: dict[str, Any] | None = None,
    slab_static_kwargs: dict[str, Any] | None = None,
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
    make_slabs_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.slabs.make_slabs_from_bulk][]
    run_static
        Whether to run the static calculation.
    slab_relax_kwargs
        Additional keyword arguments to pass to [quacc.recipes.emt.core.relax_job][].
    slab_static_kwargs
        Additional keyword arguments to pass to [quacc.recipes.emt.core.static_job][].

    Returns
    -------
    list[RunSchema | OptSchema]
        [RunSchema][quacc.schemas.ase.summarize_run] or
        [OptSchema][quacc.schemas.ase.summarize_opt_run] for each slab.
    """
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    make_slabs_kwargs = make_slabs_kwargs or {}

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    @job
    def _make_slabs(atoms):
        return make_slabs_from_bulk(atoms, **make_slabs_kwargs)

    @subflow
    def _relax_distributed(slabs):
        return [relax_job(slab, **slab_relax_kwargs) for slab in slabs]

    @subflow
    def _relax_and_static_distributed(slabs):
        return [
            static_job(
                relax_job(slab, **slab_relax_kwargs)["atoms"],
                **slab_static_kwargs,
            )
            for slab in slabs
        ]

    slabs = _make_slabs(atoms)

    return (
        _relax_and_static_distributed(slabs)
        if run_static
        else _relax_distributed(slabs)
    )
