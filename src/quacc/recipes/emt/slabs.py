"""Slab recipes for EMT."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
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
    custom_relax_job: Job | None = None,
    custom_static_job: Job | None = None,
    run_static: bool = True,
    make_slabs_kwargs: dict[str, Any] | None = None,
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
    custom_relax_job
        The relaxation job, which defaults to [quacc.recipes.emt.core.relax_job][].
    custom_static_job
        The static job, which defaults to [quacc.recipes.emt.core.static_job][].
    run_static
        Whether to run static calculations.
    make_slabs_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.slabs.make_slabs_from_bulk][]

    Returns
    -------
    list[RunSchema | OptSchema]
        [RunSchema][quacc.schemas.ase.summarize_run] or
        [OptSchema][quacc.schemas.ase.summarize_opt_run] for each slab.
    """

    return bulk_to_slabs_subflow(
        atoms,
        relax_job if custom_relax_job is None else custom_relax_job,
        static_job=(static_job if custom_static_job is None else custom_relax_job)
        if run_static
        else None,
        make_slabs_kwargs=make_slabs_kwargs,
    )
