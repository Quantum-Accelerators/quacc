"""Slab recipes for EMT"""
from __future__ import annotations

import covalent as ct
from ase import Atoms

from quacc.recipes.emt.core import relax_job, static_job
from quacc.util.slabs import make_max_slabs_from_bulk


def bulk_to_slabs_flow(
    atoms: Atoms,
    slabgen_kwargs: dict | None = None,
    slab_relax_electron: ct.electron = relax_job,
    slab_static_electron: ct.electron | None = static_job,
    slab_relax_kwargs: dict | None = None,
    slab_static_kwargs: dict | None = None,
) -> list[dict]:
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
    list[dict]
        List of dictionary of results from quacc.schemas.ase.summarize_run or quacc.schemas.ase.summarize_opt_run
    """

    slab_relax_kwargs = slab_relax_kwargs or {}
    slab_static_kwargs = slab_static_kwargs or {}
    slabgen_kwargs = slabgen_kwargs or {}

    if "relax_cell" not in slab_relax_kwargs:
        slab_relax_kwargs["relax_cell"] = False

    @ct.electron
    @ct.lattice
    def _relax_distributed(slabs):
        return [slab_relax_electron(slab, **slab_relax_kwargs) for slab in slabs]

    @ct.electron
    @ct.lattice
    def _relax_and_static_distributed(slabs):
        return [
            slab_static_electron(
                slab_relax_electron(slab, **slab_relax_kwargs)["atoms"],
                **slab_static_kwargs,
            )
            for slab in slabs
        ]

    slabs = ct.electron(make_max_slabs_from_bulk)(atoms, **slabgen_kwargs)

    if slab_static_electron is None:
        return _relax_distributed(slabs)
    else:
        return _relax_and_static_distributed(slabs)
