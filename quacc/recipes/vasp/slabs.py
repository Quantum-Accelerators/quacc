"""Recipes for slabs"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import covalent as ct
from ase.atoms import Atoms
from covalent._workflow.electron import Electron

from quacc.calculators.vasp import Vasp
from quacc.schemas.vasp import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc
from quacc.util.slabs import make_max_slabs_from_bulk, slab_to_adsorbates


def slab_static_job(
    atoms: Atoms, preset: str | None = None, swaps: dict[str, Any] | None = None
) -> dict[str, Any]:
    """
    Function to carry out a single-point calculation on a slab.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        dictionary of custom kwargs for the calculator.

    Returns
    -------
    summary
        Dictionary of the run summary.
    """

    swaps = swaps or {}

    defaults = {
        "auto_dipole": True,
        "ismear": -5,
        "laechg": True,
        "lcharg": True,
        "lvhar": True,
        "lwave": True,
        "nedos": 5001,
        "nsw": 0,
    }
    flags = merge_dicts(defaults, swaps)

    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms)
    summary = summarize_run(atoms)

    return summary


def slab_relax_job(
    atoms: Atoms, preset: str | None = None, swaps: dict[str, Any] | None = None
) -> dict[str, Any]:
    """
    Function to relax a slab.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    summary
        Dictionary of the run summary.
    """

    swaps = swaps or {}

    defaults = {
        "auto_dipole": True,
        "ediffg": -0.02,
        "isif": 2,
        "ibrion": 2,
        "isym": 0,
        "lcharg": False,
        "lwave": False,
        "nsw": 200,
    }
    flags = merge_dicts(defaults, swaps)

    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms)
    summary = summarize_run(atoms)

    return summary


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

    relax_electron: Electron | None = ct.electron(slab_relax_job)
    static_electron: Electron | None = ct.electron(slab_static_job)
    relax_kwargs: dict[str, Any] | None = None
    static_kwargs: dict[str, Any] | None = None

    def run(
        self,
        atoms: Atoms,
        max_slabs: int = None,
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

        slabs = make_max_slabs_from_bulk(atoms, max_slabs=max_slabs, **slabgen_kwargs)

        if self.relax_electron and self.static_electron:
            return _relax_and_static_distributed(slabs)
        elif self.relax_electron:
            return _relax_distributed(slabs)
        elif self.static_electron:
            return _static_distributed(slabs)


@dataclass
class SlabToAdsFlow:
    """
    Workflow consisting of:
    1. Slab-adsorbate generation
    2. Slab-adsorbate relaxations (optional)
    3. Slab-adsorbate statics (optional)

    Parameters
    ----------
    relax_electron
        Default to use for the relaxation of the slab structure.
    static_electron
        Default to use for the static calculation of the slab structures.
    relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    static_kwargs
        Additional keyword arguments to pass to the static calculation.
    """

    relax_electron: Electron | None = ct.electron(slab_relax_job)
    static_electron: Electron | None = ct.electron(slab_static_job)
    relax_kwargs: dict[str, Any] | None = None
    static_kwargs: dict[str, Any] | None = None

    def run(
        self,
        slab: Atoms,
        adsorbate: Atoms,
        **make_ads_kwargs,
    ):
        """
        Make the run.

        Parameters
        ----------
        slab
            .Atoms object for the slab structure.
        adsorbate
            .Atoms object for the adsorbate.
        **make_ads_kwargs
            Additional keyword arguments to pass to make_adsorbate_structures()
        """

        relax_kwargs = self.relax_kwargs or {}
        static_kwargs = self.static_kwargs or {}
        make_ads_kwargs = make_ads_kwargs or {}

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

        ads_slabs = slab_to_adsorbates(slab, adsorbate, **make_ads_kwargs)

        if self.relax_electron and self.static_electron:
            return _relax_and_static_distributed(ads_slabs)
        elif self.relax_electron:
            return _relax_distributed(ads_slabs)
        elif self.static_electron:
            return _static_distributed(ads_slabs)
