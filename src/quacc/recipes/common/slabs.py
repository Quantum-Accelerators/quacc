"""Common slab workflows"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import subflow
from quacc.atoms.slabs import make_slabs_from_bulk

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema


def common_bulk_to_slabs_flow(
    atoms: Atoms,
    relax_job: Callable,
    static_job: Callable | None,
    make_slabs_kwargs: dict[str, Any] | None = None,
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
    relax_job
        The relaxation function.
    static_job
        The static function.
    make_slabs_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.slabs.make_slabs_from_bulk][]
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

    @subflow
    def _relax_job_distributed(atoms: Atoms) -> list:
        slabs = make_slabs_from_bulk(atoms, **make_slabs_kwargs)
        return [relax_job(slab, **slab_relax_kwargs) for slab in slabs]

    @subflow
    def _relax_and_static_job_distributed(atoms: Atoms) -> list:
        slabs = make_slabs_from_bulk(atoms, **make_slabs_kwargs)
        return [
            static_job(
                relax_job(slab, **slab_relax_kwargs)["atoms"],
                **slab_static_kwargs,
            )
            for slab in slabs
        ]

    return (
        _relax_and_static_job_distributed(atoms)
        if static_job
        else _relax_job_distributed(atoms)
    )
