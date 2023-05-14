from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import covalent as ct
from ase.atoms import Atoms
from covalent._workflow.electron import Electron

from quacc.recipes.emt.core import relax_job, static_job
from quacc.util.slabs import make_max_slabs_from_bulk


@dataclass
class BulkToSlabsFlow:
    """
    Workflow consisting of:
    1. Slab generation
    2. Slab relaxations (optional)
    3. Slab statics (optional)

    Parameters
    ----------
    relax_electron
        Default to use for the relaxation of the slab structures.
    static_electron
        Default to use for the static calculation of the slab structures.
    relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    static_kwargs
        Additional keyword arguments to pass to the static calculation.
    """

    relax_electron: Electron | None = ct.electron(relax_job)
    static_electron: Electron | None = ct.electron(static_job)
    relax_kwargs: dict[str, Any] | None = None
    static_kwargs: dict[str, Any] | None = None

    def run(
        self,
        atoms: Atoms,
        slabgen_kwargs: dict[str, Any] = None,
    ):
        """
        Make the workflow.

        Parameters
        ----------
        atoms
            .Atoms object for the structure.
        slabgen_kwargs
            Additional keyword arguments to pass to make_max_slabs_from_bulk()
        """

        relax_kwargs = self.relax_kwargs or {}
        static_kwargs = self.static_kwargs or {}
        slabgen_kwargs = slabgen_kwargs or {}

        if not self.relax_electron and not self.static_electron:
            raise ValueError(
                "At least one of relax_electron or static_electron must be defined."
            )

        @ct.electron
        @ct.lattice
        def _relax_distributed(slabs):
            return [self.relax_electron(slab, **relax_kwargs) for slab in slabs]

        @ct.electron
        @ct.lattice
        def _static_distributed(slabs):
            return [self.static_electron(slab, **static_kwargs) for slab in slabs]

        @ct.electron
        @ct.lattice
        def _relax_and_static_distributed(slabs):
            return [
                self.static_electron(
                    self.relax_electron(slab, **relax_kwargs)["atoms"], **static_kwargs
                )
                for slab in slabs
            ]

        slabs = make_max_slabs_from_bulk(atoms, **slabgen_kwargs)

        if self.relax_electron and self.static_electron:
            return _relax_and_static_distributed(slabs)
        elif self.relax_electron:
            return _relax_distributed(slabs)
        elif self.static_electron:
            return _static_distributed(slabs)
