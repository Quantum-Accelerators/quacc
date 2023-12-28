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
    decorators: dict[str, Callable | None] | None = None,
    parameters: dict[str, Any] | None = None,
) -> list[RunSchema | OptSchema]:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Optional slab statics

    Parameters
    ----------
    atoms
        Atoms object
    make_slabs_kwargs
        Additional keyword arguments to pass to [quacc.atoms.slabs.make_slabs_from_bulk][]
    decorators
        Custom decorators to apply to each Job in the Flow.
    parameters
        Custom parameters to pass to each Job in the Flow.

    Returns
    -------
    list[RunSchema | OptSchema]
        [RunSchema][quacc.schemas.ase.summarize_run] or
        [OptSchema][quacc.schemas.ase.summarize_opt_run] for each slab.
    """
    slab_relax_job_, slab_static_job_, bulk_to_slabs_subflow_ = customize_funcs(
        {
            "slab_relax_job": relax_job,
            "slab_static_job": static_job,
            "bulk_to_slabs_subflow": bulk_to_slabs_subflow,
        },
        decorators=decorators,
        parameters=parameters,
    )

    return bulk_to_slabs_subflow_(
        atoms,
        slab_relax_job_,
        static_job=slab_static_job_,
        make_slabs_kwargs=make_slabs_kwargs,
    )
