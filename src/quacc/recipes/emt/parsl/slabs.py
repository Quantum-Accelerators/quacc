"""Slab recipes for EMT"""
from __future__ import annotations

from ase import Atoms
from parsl import python_app
from parsl.app.python import PythonApp
from parsl.dataflow.futures import AppFuture

from quacc.recipes.emt.core import relax_job, static_job
from quacc.schemas.ase import OptSchema, RunSchema
from quacc.util.slabs import make_max_slabs_from_bulk


# See https://github.com/Parsl/parsl/issues/2793 for why we need to strip the @ct.electron
# decorator off the PythonApp kwargs
def bulk_to_slabs_flow(
    atoms: Atoms | dict,
    make_slabs_kwargs: dict | None = None,
    slab_relax: PythonApp = python_app(relax_job.electron_object.function),
    slab_static: PythonApp | None = python_app(static_job.electron_object.function),
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> AppFuture[RunSchema | OptSchema]:
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
    AppFuture[RunSchema | OptSchema]
        List of AppFuture objects, each of which resolves to a dictionary of results
        from quacc.schemas.ase.summarize_run or quacc.schemas.ase.summarize_opt_run
    """
    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    make_slabs_kwargs = make_slabs_kwargs or {}

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    @python_app
    def _make_slabs(atoms):
        atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
        return make_max_slabs_from_bulk(atoms, **make_slabs_kwargs)

    slabs = _make_slabs(atoms).result()

    futures = []
    for slab in slabs:
        slab_relax_future = slab_relax(slab, **slab_relax_kwargs)

        if slab_static:
            slab_static_future = slab_static(slab_relax_future, **slab_static_kwargs)
            futures.append(slab_static_future)
        else:
            futures.append(slab_relax_future)

    return futures
