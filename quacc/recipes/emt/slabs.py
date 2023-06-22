"""Slab recipes for EMT"""
from __future__ import annotations

from dataclasses import dataclass

import covalent as ct
from ase import Atoms
from covalent._workflow.electron import Electron

from quacc.recipes.emt.core import relax_job, static_job
from quacc.util.slabs import make_max_slabs_from_bulk


@dataclass
class BulkToSlabsFlow:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations

    3. Slab statics (optional)

    Parameters
    ----------
    slab_relax_electron
        Default Electron to use for the relaxation of the slab structures.
    slab_static_electron
        Default Electron to use for the static calculation of the slab structures.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.
    """

    slab_relax_electron: Electron = relax_job
    slab_static_electron: Electron | None = static_job
    slab_relax_kwargs: dict | None = None
    slab_static_kwargs: dict | None = None

    def run(
        self,
        atoms: Atoms,
        slabgen_kwargs: dict | None = None,
    ) -> list[dict]:
        """
        Make the workflow.

        Parameters
        ----------
        atoms
            Atoms object for the structure.
        slabgen_kwargs
            Additional keyword arguments to pass to make_max_slabs_from_bulk()

        Returns
        -------
        list[dict]
            List of dictionary of results from quacc.schemas.ase.summarize_run or quacc.schemas.ase.summarize_opt_run
        """

        self.slab_relax_kwargs = self.slab_relax_kwargs or {}
        self.slab_static_kwargs = self.slab_static_kwargs or {}
        slabgen_kwargs = slabgen_kwargs or {}

        if "relax_cell" not in self.slab_relax_kwargs:
            self.slab_relax_kwargs["relax_cell"] = False

        @ct.electron
        @ct.lattice
        def _relax_distributed(slabs):
            return [
                self.slab_relax_electron(slab, **self.slab_relax_kwargs)
                for slab in slabs
            ]

        @ct.electron
        @ct.lattice
        def _relax_and_static_distributed(slabs):
            return [
                self.slab_static_electron(
                    self.slab_relax_electron(slab, **self.slab_relax_kwargs)["atoms"],
                    **self.slab_static_kwargs,
                )
                for slab in slabs
            ]

        slabs = ct.electron(make_max_slabs_from_bulk)(atoms, **slabgen_kwargs)

        if self.slab_static_electron is None:
            return _relax_distributed(slabs)
        else:
            return _relax_and_static_distributed(slabs)
