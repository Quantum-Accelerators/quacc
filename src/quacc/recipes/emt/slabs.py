"""Slab recipes for EMT."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.common.slabs import common_bulk_to_slabs_flow
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

    slab_relax_kwargs = slab_relax_kwargs or {}
    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    return common_bulk_to_slabs_flow(
        atoms,
        relax_job,
        static_job if run_static else None,
        make_slabs_kwargs=make_slabs_kwargs,
        slab_relax_kwargs=slab_relax_kwargs,
        slab_static_kwargs=slab_static_kwargs,
    )
