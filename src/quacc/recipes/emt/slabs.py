"""Slab recipes for EMT."""
from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

from quacc import flow
from quacc.atoms.slabs import make_slabs_from_bulk
from quacc.recipes.common.slabs import bulk_to_slabs_subflow
from quacc.recipes.emt.core import relax_job, static_job

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc import Job
    from quacc.schemas._aliases.ase import OptSchema, RunSchema


@flow
def bulk_to_slabs_flow(
    atoms: Atoms,
    make_slabs_kwargs: dict[str, Any] | None = None,
    slab_relax_job: Job | None = relax_job,
    slab_static_job: Job | None = static_job,
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

    make_slabs_kwargs = make_slabs_kwargs or {}

    return bulk_to_slabs_subflow(
        atoms,
        slab_relax_job,
        static_job=static_job if slab_static_job else None,
        make_slabs_fn=partial(make_slabs_from_bulk, **make_slabs_kwargs),
    )
