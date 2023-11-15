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

    from ase import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema


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

    relax_fn = partial(relax_job, **slab_relax_kwargs if slab_relax_kwargs else {})
    static_fn = (
        partial(static_job, **slab_static_kwargs if slab_static_kwargs else {})
        if run_static
        else None
    )
    make_slabs_fn = partial(
        make_slabs_from_bulk, **make_slabs_kwargs if make_slabs_kwargs else {}
    )

    return bulk_to_slabs_subflow(
        atoms, relax_fn, static_fn=static_fn, make_slabs_fn=make_slabs_fn
    )
