"""Recipes for slabs"""
from __future__ import annotations

from dataclasses import dataclass

import covalent as ct
from ase import Atoms
from covalent._workflow.electron import Electron

from quacc.calculators.vasp import Vasp
from quacc.schemas.vasp import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import merge_dicts
from quacc.util.slabs import make_adsorbate_structures, make_max_slabs_from_bulk


@ct.electron
def slab_static_job(atoms: Atoms, preset: str = None, swaps: dict = None) -> dict:
    """
    Function to carry out a single-point calculation on a slab.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use.
    swaps
        dictionary of custom kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.vasp.summarize_run
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

    return summarize_run(atoms, additional_fields={"name": "VASP Slab Static"})


@ct.electron
def slab_relax_job(atoms: Atoms, preset: str = None, swaps: dict = None) -> dict:
    """
    Function to relax a slab.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.vasp.summarize_run
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

    return summarize_run(atoms, additional_fields={"name": "VASP Slab Relax"})


@dataclass
class BulkToSlabsFlow:
    """
    Workflow consisting of:

    1. Slab generation

    2. Slab relaxations (optional)

    3. Slab statics (optional)

    Parameters
    ----------
    slab_relax_electron
        Default to use for the relaxation of the slab structures.
    slab_static_electron
        Default to use for the static calculation of the slab structures.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.
    """

    slab_relax_electron: Electron | None = slab_relax_job
    slab_static_electron: Electron | None = slab_static_job
    slab_relax_kwargs: dict = None
    slab_static_kwargs: dict = None

    def run(
        self,
        atoms: Atoms,
        slabgen_kwargs: dict = None,
    ) -> list[dict]:
        """
        Run the workflow.

        Parameters
        ----------
        atoms
            Atoms object for the structure.
        slabgen_kwargs
            Additional keyword arguments to pass to make_max_slabs_from_bulk()

        Returns
        -------
        list[dict]
            List of dictionary results from quacc.schemas.vasp.summarize_run
        """

        self.relax_kwargs = self.slab_relax_kwargs or {}
        self.static_kwargs = self.slab_static_kwargs or {}
        slabgen_kwargs = slabgen_kwargs or {}

        if not self.slab_relax_electron and not self.slab_static_electron:
            raise ValueError(
                "At least one of slab_relax_electron or slab_static_electron must be defined."
            )

        @ct.electron
        @ct.lattice
        def _relax_distributed(slabs):
            return [
                self.slab_relax_electron(slab, **self.relax_kwargs) for slab in slabs
            ]

        @ct.electron
        @ct.lattice
        def _static_distributed(slabs):
            return [
                self.slab_static_electron(slab, **self.static_kwargs) for slab in slabs
            ]

        @ct.electron
        @ct.lattice
        def _relax_and_static_distributed(slabs):
            return [
                self.slab_static_electron(
                    self.slab_relax_electron(slab, **self.relax_kwargs)["atoms"],
                    **self.static_kwargs,
                )
                for slab in slabs
            ]

        slabs = ct.electron(make_max_slabs_from_bulk)(atoms, **slabgen_kwargs)

        if self.slab_relax_electron and self.slab_static_electron:
            return _relax_and_static_distributed(slabs)
        elif self.slab_relax_electron:
            return _relax_distributed(slabs)
        elif self.slab_static_electron:
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
    slab_relax_electron
        Default to use for the relaxation of the slab structure.
    slab_static_electron
        Default to use for the static calculation of the slab structures.
    slab_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    slab_static_kwargs
        Additional keyword arguments to pass to the static calculation.
    """

    slab_relax_electron: Electron | None = ct.electron(slab_relax_job)
    slab_static_electron: Electron | None = ct.electron(slab_static_job)
    slab_relax_kwargs: dict = None
    slab_static_kwargs: dict = None

    def run(
        self,
        slab: Atoms,
        adsorbate: Atoms,
        **make_ads_kwargs,
    ) -> dict:
        """
        Make the run.

        Parameters
        ----------
        slab
            Atoms object for the slab structure.
        adsorbate
            Atoms object for the adsorbate.
        **make_ads_kwargs
            Additional keyword arguments to pass to make_adsorbate_structures()
        """

        self.slab_relax_kwargs = self.slab_relax_kwargs or {}
        self.slab_static_kwargs = self.slab_static_kwargs or {}
        make_ads_kwargs = make_ads_kwargs or {}

        if not self.slab_relax_electron and not self.slab_static_electron:
            raise ValueError(
                "At least one of slab_relax_electron or slab_static_electron must be defined."
            )

        @ct.electron
        @ct.lattice
        def _relax_distributed(slabs):
            return [
                self.slab_relax_electron(slab, **self.slab_relax_kwargs)
                for slab in slabs
            ]

        @ct.electron
        @ct.lattice
        def _static_distributed(slabs):
            return [
                self.slab_static_electron(slab, **self.slab_static_kwargs)
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

        ads_slabs = ct.electron(make_adsorbate_structures)(
            slab, adsorbate, **make_ads_kwargs
        )

        if self.slab_relax_electron and self.slab_static_electron:
            return _relax_and_static_distributed(ads_slabs)
        elif self.slab_relax_electron:
            return _relax_distributed(ads_slabs)
        elif self.slab_static_electron:
            return _static_distributed(ads_slabs)
