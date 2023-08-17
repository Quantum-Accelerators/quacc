"""Slab recipes for EMT"""
from __future__ import annotations

from typing import TYPE_CHECKING

from parsl import join_app, python_app

from quacc.recipes.emt.core import relax_job, static_job
from quacc.schemas.atoms import fetch_atoms
from quacc.util.slabs import make_max_slabs_from_bulk

if TYPE_CHECKING:
    from ase import Atoms
    from parsl.app.python import PythonApp
    from parsl.dataflow.futures import AppFuture

    from quacc.schemas.ase import OptSchema, RunSchema


def bulk_to_slabs_flow(
    atoms: Atoms | dict,
    make_slabs_kwargs: dict | None = None,
    slab_relax: PythonApp = python_app(relax_job.electron_object.function),
    slab_static: PythonApp | None = python_app(static_job.electron_object.function),
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> AppFuture[list[RunSchema | OptSchema]]:
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
        Additional keyword arguments to pass to make_max_slabs_from_bulk()
    slab_relax
        Default PythonApp to use for the relaxation of the slab structures.
    slab_static
        Default PythonApp to use for the static calculation of the slab structures.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    AppFuture[list[RunSchema | OptSchema]]
        An AppFuture whose .result() is a list of dictionary of results from
        quacc.schemas.ase.summarize_run or quacc.schemas.ase.summarize_opt_run
    """
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    make_slabs_kwargs = make_slabs_kwargs or {}

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    @python_app
    def _make_slabs(atoms):
        atoms = fetch_atoms(atoms)
        return make_max_slabs_from_bulk(atoms, **make_slabs_kwargs)

    @join_app
    def _relax_distributed(slabs):
        return [slab_relax(slab, **slab_relax_kwargs) for slab in slabs]

    @join_app
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
