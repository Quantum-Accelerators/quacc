"""Slab recipes for EMT"""
from __future__ import annotations

from ase import Atoms
from parsl import join_app, python_app
from parsl.app.python import PythonApp
from parsl.dataflow.futures import AppFuture

from quacc.recipes.emt.core import relax_job, static_job


@join_app
def bulk_to_slabs_app(
    atoms: Atoms,
    slabgen_kwargs: dict | None = None,
    slab_relax_app: PythonApp = python_app(relax_job),
    slab_static_app: PythonApp | None = python_app(static_job),
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> AppFuture:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    slabgen_kwargs
        Additional keyword arguments to pass to make_max_slabs_from_bulk()
    slab_relax_electron
        Default Electron to use for the relaxation of the slab structures.
    slab_static_electron
        Default Electron to use for the static calculation of the slab structures.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    AppFuture
        An AppFuture whose .result() is a list[dict]
    """

    from quacc.util.slabs import make_max_slabs_from_bulk

    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    slabgen_kwargs = slabgen_kwargs or {}

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    def _relax_distributed(slabs):
        return [slab_relax_app(slab, **slab_relax_kwargs) for slab in slabs]

    def _relax_and_static_distributed(slabs):
        return [
            slab_static_app(
                slab_relax_app(slab, **slab_relax_kwargs).result()["atoms"],
                **slab_static_kwargs,
            )
            for slab in slabs
        ]

    slabs = make_max_slabs_from_bulk(atoms, **slabgen_kwargs)

    if slab_relax_app and slab_static_app:
        return _relax_and_static_distributed(slabs)
    else:
        return _relax_distributed(slabs)
