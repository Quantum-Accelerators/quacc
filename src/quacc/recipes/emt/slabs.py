"""Slab recipes for EMT."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.common.slabs import bulk_to_slabs_subflow
from quacc.recipes.emt.core import relax_job, static_job
from quacc.wflow_tools.decorators import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc import Job
    from quacc.schemas._aliases.ase import OptSchema, RunSchema


@flow
def bulk_to_slabs_flow(
    atoms: Atoms,
    make_slabs_kwargs: dict[str, Any] | None = None,
    slab_relax_job: Job = relax_job,
    slab_static_job: Job | None = static_job,
    decorators: dict[str, Callable | None] | None = None,
    common_kwargs: dict[str, Any] | None = None,
) -> list[RunSchema | OptSchema]:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Optional Slab statics

    Parameters
    ----------
    atoms
        Atoms object
    make_slabs_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.slabs.make_slabs_from_bulk][]
    slab_relax_job
        The slab relaxation job, which defaults to [quacc.recipes.emt.core.relax_job][].
    slab_static_job
        The slab static job, which defaults to [quacc.recipes.emt.core.static_job][].
    decorators
        Custom decorators to apply to each Job in the Flow.
    common_kwargs
        Common keyword arguments to pass to each Job in the Flow.

    Returns
    -------
    list[RunSchema | OptSchema]
        [RunSchema][quacc.schemas.ase.summarize_run] or
        [OptSchema][quacc.schemas.ase.summarize_opt_run] for each slab.
    """
    slab_relax_job_, slab_static_job_, bulk_to_slabs_subflow_ = customize_funcs(
        {
            "slab_relax_job": slab_relax_job,
            "slab_static_job": slab_static_job,
            "bulk_to_slabs_subflow": bulk_to_slabs_subflow,
        },
        decorators=decorators,
        common_kwargs=common_kwargs,
    )

    return bulk_to_slabs_subflow_(
        atoms,
        slab_relax_job_,
        static_job=slab_static_job_,
        make_slabs_kwargs=make_slabs_kwargs,
    )
