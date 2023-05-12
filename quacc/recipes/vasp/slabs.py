"""Recipes for slabs"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any, List

import covalent as ct
from ase.atoms import Atoms
from covalent._workflow.electron import Electron

from quacc.calculators.vasp import Vasp
from quacc.recipes.vasp.core import relax_job, static_job
from quacc.schemas.vasp import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc
from quacc.util.slabs import (
    get_surface_energy,
    make_adsorbate_structures,
    make_max_slabs_from_bulk,
    slab_to_adsorbates,
)


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
    """

    relax_electron: Electron | None = ct.electron(slab_relax_job)
    static_electron: Electron | None = ct.electron(slab_static_job)

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

        slabgen_kwargs = slabgen_kwargs or {}

        @ct.electron
        @ct.lattice
        def _relax_distributed(slabs):
            return [self.relax_electron(slab) for slab in slabs]

        @ct.electron
        @ct.lattice
        def _relax_and_static_distributed(slabs):
            return [
                self.static_electron(self.relax_electron(slab)["atoms"])
                for slab in slabs
            ]

        slabs = make_max_slabs_from_bulk(atoms, max_slabs=max_slabs, **slabgen_kwargs)

        if self.relax_electron and self.static_electron:
            return _relax_and_static_distributed(slabs)
        elif self.relax_electron:
            return _relax_distributed(slabs)
        else:
            return slabs


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
        Default to use for the relaxation of the slab structures.
    static_electron
        Default to use for the static calculation of the slab structures.
    """

    relax_electron: Electron | None = ct.electron(slab_relax_job)
    static_electron: Electron | None = ct.electron(slab_static_job)

    def run(self, slab: Atoms, adsorbate: Atoms, **make_ads_kwargs):
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

        make_ads_kwargs = make_ads_kwargs or {}

        @ct.electron
        @ct.lattice
        def _relax_distributed(slabs):
            return [self.relax_electron(slab) for slab in slabs]

        @ct.electron
        @ct.lattice
        def _relax_and_static_distributed(slabs):
            return [
                self.static_electron(self.relax_electron(slab)["atoms"])
                for slab in slabs
            ]

        ads_slabs = make_adsorbate_structures(slab, adsorbate, **make_ads_kwargs)

        if self.relax_electron and self.static_electron:
            return _relax_and_static_distributed(ads_slabs)
        elif self.relax_electron:
            return _relax_distributed(ads_slabs)
        else:
            return ads_slabs


def _get_slab_stability(
    bulk_summary: dict[str, Any],
    slab_summaries: dict[str, Any],
    n_stable_slabs: int = 1,
) -> dict[str, Any]:
    """
    A function that determine the most stable surface slabs (based on surface energy) for
    a given bulk summary and list of slab summaries.

    Parameters
    ----------
    bulk_summary
        Output of a VASP job corresponding to the bulk structure.
    slab_summaries
        list of outputs of VASP jobs corresponding to the slab structures.
    n_stable_slabs
        The n most stable slabs are returned.

    Returns
    -------
    dict
        VASP output summaries for the stable and unstable slabs formatted as
        {
            {"stable_slabs": {"all_atoms": [.Atoms, .Atoms, ...], "all_outputs": [...],
            {"unstable_slabs": {"all_atoms": [.Atoms, .Atoms, ...], "all_outputs": [...]
        }
    """
    bulk = bulk_summary["atoms"]
    bulk_energy = bulk_summary["output"]["energy"]
    surface_energies = []

    if n_stable_slabs > len(slab_summaries):
        warnings.warn(
            "n_stable_slabs is larger than the number of slabs. Setting n_stable_slabs to the number of slabs."
        )
        n_stable_slabs = len(slab_summaries)

    # Iterate through each slab summary and determine the most stable slab
    for slab_summary in slab_summaries:
        slab = slab_summary["atoms"]
        slab_energy = slab_summary["output"]["energy"]

        # Calculate the surface energy
        surface_energy = get_surface_energy(bulk, slab, bulk_energy, slab_energy)

        # Insert the surface energy into the slab summary
        slab_summary["avg_surface_energy"] = surface_energy

        # Store the slab energy in a convenient dict
        surface_energies.append(surface_energy)

    slab_summaries_sorted = [
        slab_summary
        for _, slab_summary in sorted(zip(surface_energies, slab_summaries))
    ]
    stable_slab_summaries = slab_summaries_sorted[0:n_stable_slabs]
    unstable_slab_summaries = slab_summaries_sorted[n_stable_slabs:]

    output = {
        "stable_slabs": {
            "all_atoms": [summary["atoms"] for summary in stable_slab_summaries],
            "all_outputs": stable_slab_summaries,
        },
        "unstable_slabs": {
            "all_atoms": [summary["atoms"] for summary in unstable_slab_summaries]
            if unstable_slab_summaries
            else [],
            "all_outputs": unstable_slab_summaries,
        },
    }

    return output
