from dataclasses import dataclass
from typing import Any, Dict

import numpy as np
from ase.atoms import Atoms
from jobflow import Flow, Maker, Response, job
from jobflow.core.flow import JobOrder

from quacc.calculators.vasp import SmartVasp
from quacc.recipes.vasp.core import DoubleRelaxJob, StaticJob
from quacc.schemas.vasp import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc
from quacc.util.slabs import (
    get_cleavage_energy,
    make_adsorbate_structures,
    make_max_slabs_from_bulk,
)


@dataclass
class SlabStaticJob(Maker):
    """
    Class to carry out a single-point calculation on a slab.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "VASP-SlabStatic"
    preset: str = None
    swaps: Dict[str, Any] = None

    @job
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object

        Returns
        -------
        Dict
            Summary of the run.
        """
        swaps = self.swaps or {}
        defaults = {
            "auto_dipole": True,
            "ismear": -5,
            "isym": 2,
            "laechg": True,
            "lcharg": True,
            "lvhar": True,
            "lwave": True,
            "nedos": 5001,
            "nsw": 0,
            "sigma": 0.05,
        }
        flags = merge_dicts(defaults, swaps, remove_none=True)

        atoms = SmartVasp(atoms, preset=self.preset, **flags)
        atoms = run_calc(atoms)
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary


@dataclass
class SlabRelaxJob(Maker):
    """
    Class to relax a slab.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "VASP-SlabRelax"
    preset: str = None
    swaps: Dict[str, Any] = None

    @job
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object

        Returns
        -------
        Dict:
            Summary of the run.
        """
        swaps = self.swaps or {}
        defaults = {
            "auto_dipole": True,
            "ediffg": -0.02,
            "isif": 2,
            "ibrion": 2,
            "ismear": 0,
            "isym": 0,
            "lcharg": False,
            "lwave": False,
            "nsw": 200,
            "sigma": 0.05,
        }
        flags = merge_dicts(defaults, swaps, remove_none=True)

        atoms = SmartVasp(atoms, preset=self.preset, **flags)
        atoms = run_calc(atoms)
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary


@dataclass
class BulkToSlabsJob(Maker):
    """
    Class to convert a bulk structure to slabs,
    along with the relaxations and statics for the slabs.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use. Applies to all jobs in the flow.
    slab_relax_job
        Maker to use for the SlabRelax job.
    slab_static_job
        Default to use for the SlabStatic job.
    swaps
        Dictionary of custom kwargs for the calculator.
        Applies to all jobs in the flow.
    """

    name: str = "VASP-BulkToSlab"
    preset: str = None
    slab_relax_job: Maker | None = SlabRelaxJob()
    slab_static_job: Maker | None = SlabStaticJob()
    swaps: Dict[str, Any] = None

    @job
    def make(self, atoms: Atoms, max_slabs: int = None, **slabgen_kwargs) -> Response:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object
        max_slabs
            Maximum number of slabs to make. None implies no upper limit.
        slabgen_kwargs
            Additional keyword arguments to pass to make_max_slabs_from_bulk()

        Returns
        -------
        Response
            A Flow of relaxation and static jobs for the generated slabs.
        """
        slabgen_kwargs = slabgen_kwargs or {}
        if self.preset:
            self.slab_static_job.preset = self.preset
            self.slab_relax_job.preset = self.preset
        if self.swaps:
            self.slab_static_job.swaps = self.swaps
            self.slab_relax_job.swaps = self.swaps

        slabs = make_max_slabs_from_bulk(atoms, max_slabs=max_slabs, **slabgen_kwargs)
        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = self.slab_relax_job.make(slab)
            jobs.append(relax_job)

            static_job = self.slab_static_job.make(relax_job.output["atoms"])
            jobs.append(static_job)
            outputs.append(static_job.output)

        return Response(replace=Flow(jobs, output=outputs, order=JobOrder.LINEAR))


@dataclass
class SlabToAdsorbatesJob(Maker):
    """
    Class to convert a slab structure to one with an adsorbate present,
    along with the relaxations and statics for the slab-adsorbate systems.
    Multiple slab-adsorbate systems will be generated, one for each unique
    binding site.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use. Applies to all jobs in the flow.
    slab_relax_job
        Maker to use for the SlabRelax job.
    slab_static_job
        Default to use for the SlabStatic job.
    swaps
        Dictionary of custom kwargs for the calculator.
        Applies to all jobs in the flow.
    """

    name: str = "VASP-SlabToAdsorbates"
    preset: str = None
    swaps: Dict[str, Any] = None
    slab_relax_job: Maker | None = SlabRelaxJob()
    slab_static_job: Maker | None = SlabStaticJob()

    @job
    def make(self, atoms: Atoms, adsorbate: Atoms, **make_ads_kwargs) -> Response:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object for the structure.
        adsorbate
            .Atoms object for the adsorbate.
        make_ads_kwargs
            Additional keyword arguments to pass to make_adsorbate_structures()

        Returns
        -------
        Response
            A Flow of relaxation and static jobs for the generated slabs with adsorbates.
        """
        make_ads_kwargs = make_ads_kwargs or {}
        if self.preset:
            self.slab_static_job.preset = self.preset
            self.slab_relax_job.preset = self.preset
        if self.swaps:
            self.slab_static_job.swaps = self.swaps
            self.slab_relax_job.swaps = self.swaps

        slabs = make_adsorbate_structures(atoms, adsorbate, **make_ads_kwargs)

        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = self.slab_relax_job.make(slab)
            jobs.append(relax_job)

            static_job = self.slab_static_job.make(relax_job.output["atoms"])
            jobs.append(static_job)
            outputs.append(static_job.output)

        return Response(replace=Flow(jobs, output=outputs, order=JobOrder.LINEAR))


@dataclass
class BulkToAdsorbatesFlow(Maker):
    """
    Maker to get adsorption energies from a bulk structure.
    """

    name: str = "VASP-BulkToAdsorbates"
    preset: str = None
    bulk_relax_job: Maker | None = DoubleRelaxJob()
    bulk_static_job: Maker | None = StaticJob()
    bulk_to_slabs_job: Maker = BulkToSlabsJob()
    slab_to_adsorbates_job: Maker = SlabToAdsorbatesJob()
    swaps: Dict[str, Any] = None

    def make(
        self,
        atoms: Atoms,
        adsorbate: Atoms,
        max_slabs: int = None,
        slabgen_kwargs: Dict[str, Any] = None,
        make_ads_kwargs: Dict[str, Any] = None,
    ) -> Flow:
        jobs = []
        slabgen_kwargs = slabgen_kwargs or {}
        make_ads_kwargs = make_ads_kwargs or {}

        if self.bulk_relax_job:
            if self.preset:
                self.bulk_relax_job.preset = self.preset
            if self.swaps:
                self.bulk_relax_job.swaps = self.swaps
            bulk_relax_job = self.bulk_relax_job.make(atoms)
            atoms = bulk_relax_job.output["atoms"]
            jobs.append(bulk_relax_job)

        if self.bulk_static_job:
            if self.preset:
                self.bulk_static_job.preset = self.preset
            if self.swaps:
                self.bulk_static_job.swaps = self.swaps
            bulk_static_job = self.bulk_static_job.make(atoms)
            atoms = bulk_static_job.output["atoms"]
            jobs.append(bulk_static_job)

        bulk_to_slabs_job = self.bulk_to_slabs_job.make(
            atoms, max_slabs=max_slabs, **slabgen_kwargs
        )
        find_stable_slab_job = _get_stable_slab_summary(
            bulk_static_job.output, bulk_to_slabs_job.output
        )
        if self.preset:
            self.slab_to_adsorbates_job.preset = self.preset
        if self.swaps:
            self.slab_to_adsorbates_job.swaps = self.swaps
        slab_to_adsorbates_job = self.slab_to_adsorbates_job.make(
            find_stable_slab_job.output["stable_slab"]["atoms"],
            adsorbate,
            **make_ads_kwargs
        )

        jobs += [bulk_to_slabs_job, find_stable_slab_job, slab_to_adsorbates_job]

        return Flow(jobs)


@job
def _get_stable_slab_summary(
    bulk_summary: Dict[str, Any], slab_summaries: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Determine the most stable surface slab (based on cleavage energy) for
    a given bulk summary and list of slab summaries.

    Parameters
    ----------
    bulk_summary
        Output of a VASP job corresponding to the bulk structure.
    slab_summaries
        List of outputs of VASP jobs corresponding to the slab structures.

    Returns
    -------
    Dict
        VASP output summaries for the stable and unstable slabs formatted as
        {"stable_slab": {"atoms": ..., "output": ...}, "unstable_slabs": [...]}
    """
    min_cleave_energy = np.inf
    bulk = bulk_summary["atoms"]
    bulk_energy = bulk_summary["output"]["energy"]
    stable_slab_idx = None

    # Iterate through each slab summary and determine the most stable slab
    for i, slab_summary in enumerate(slab_summaries):
        slab = slab_summary["atoms"]
        slab_energy = slab_summary["output"]["energy"]

        # Calculate the cleavage energy
        cleave_energy = get_cleavage_energy(bulk, slab, bulk_energy, slab_energy)

        # Insert the cleavage energy into the slab summary
        slab_summary["cleavage_energy"] = cleave_energy

        # Determine if we are at a more stable slab
        if cleave_energy < min_cleave_energy:
            min_cleave_energy = cleave_energy
            stable_slab_idx = i

    # Here, we return the summary dictionaries for the stable and unstable slabs
    # with the cleavage energy inserted into the summary.
    output = {
        "stable_slab": slab_summaries[stable_slab_idx],
        "unstable_slabs": [
            summary for i, summary in enumerate(slab_summaries) if i != stable_slab_idx
        ],
    }

    return output
