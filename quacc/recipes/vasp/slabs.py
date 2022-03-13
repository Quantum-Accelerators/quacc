from dataclasses import dataclass
from typing import Any, Dict, List

import numpy as np
from ase.atoms import Atoms
from jobflow import Flow, Maker, Response, job
from jobflow.core.flow import JobOrder

from quacc.calculators.vasp import SmartVasp
from quacc.recipes.vasp.core import RelaxJob, StaticJob
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

    name: str = "VASP-BulkToSlabs"
    preset: str = None
    slab_relax_job: Maker = SlabRelaxJob()
    slab_static_job: Maker = SlabStaticJob()
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
            self.slab_relax_job.preset = self.preset
            self.slab_static_job.preset = self.preset
        if self.swaps:
            self.slab_relax_job.swaps = self.swaps
            self.slab_static_job.swaps = self.swaps

        # Generate all the slab
        slabs = make_max_slabs_from_bulk(atoms, max_slabs=max_slabs, **slabgen_kwargs)

        # Generate the jobs for each slab
        jobs = []
        outputs = []
        all_atoms = []
        for slab in slabs:
            relax_job = self.slab_relax_job.make(slab)
            static_job = self.slab_static_job.make(relax_job.output["atoms"])

            jobs += [relax_job, static_job]
            outputs.append(static_job.output)
            all_atoms.append(static_job.output["atoms"])

        return Response(
            replace=Flow(
                jobs,
                output={"all_atoms": all_atoms, "all_outputs": outputs},
                name=self.name,
                order=JobOrder.LINEAR,
            )
        )


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
    slab_relax_job: Maker = SlabRelaxJob()
    slab_static_job: Maker = SlabStaticJob()

    @job
    def make(
        self, atoms: Atoms | List[Atoms], adsorbate: Atoms, **make_ads_kwargs
    ) -> Response:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object for the slab structure. Also takes a list of Atoms objects
            for the creation of a series of slabs with adsorbates.
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

        if isinstance(atoms, Atoms):
            atoms_list = [atoms]
        else:
            atoms_list = atoms

        jobs = []
        outputs = []
        all_atoms = []
        for atoms in atoms_list:

            # Make slab-adsorbate systems
            slabs = make_adsorbate_structures(atoms, adsorbate, **make_ads_kwargs)

            # Make a relaxation+static job for each slab-adsorbate ysstem
            for slab in slabs:
                relax_job = self.slab_relax_job.make(slab)
                static_job = self.slab_static_job.make(relax_job.output["atoms"])

                jobs += [relax_job, static_job]
                outputs.append(static_job.output)
                all_atoms.append(static_job.output["atoms"])

        return Response(
            replace=Flow(
                jobs,
                output={"all_atoms": all_atoms, "all_outputs": outputs},
                name=self.name,
                order=JobOrder.LINEAR,
            )
        )


@dataclass
class BulkToAdsorbatesFlow(Maker):
    """
    Flow consisting of:
    1. Bulk relaxation (optional)
    2. Bulk static (optional)
    3. Slab generation
    4. Selection of the most stable slab (optional)
    5. Addition of adsorbates to the slabs or most stable slab
    6. Slab relaxation(s) with adsorbates at multiple binding sites
    7. Slab static(s) with adsorbates at multiple binding sites

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use. Applies to all jobs in the flow.
    bulk_relax_job
        Maker to use for the RelaxJob.
    bulk_static_job
        Default to use for the StaticJob.
    bulk_to_slabs_job
        Maker to use for the BulkToSlabsJob.
    slab_to_adsorbate_job
        Maker to use for the SlabToAdsorbatesJob.
    swaps
        Dictionary of custom kwargs for the calculator.
        Applies to all jobs in the flow.
    """

    name: str = "VASP-BulkToAdsorbates"
    preset: str = None
    bulk_relax_job: Maker | None = RelaxJob()
    bulk_static_job: Maker | None = StaticJob()
    bulk_to_slabs_job: Maker = BulkToSlabsJob()
    slab_to_adsorbates_job: Maker = SlabToAdsorbatesJob()
    swaps: Dict[str, Any] = None

    def make(
        self,
        atoms: Atoms,
        adsorbate: Atoms,
        stable_slab: bool = True,
        max_slabs: int = None,
        slabgen_kwargs: Dict[str, Any] = None,
        make_ads_kwargs: Dict[str, Any] = None,
    ) -> Flow:
        """
        Make the Flow.

        Parameters
        ----------
        atoms
            .Atoms object for the structure.
        adsorbate
            .Atoms object for the adsorbate.
        stable_slab
            Whether to add adsorbates to only the most stable slab.
        max_slabs
            Maximum number of slabs to make. None implies no upper limit.
        slabgen_kwargs
            Additional keyword arguments to pass to make_max_slabs_from_bulk()
        make_ads_kwargs
            Additional keyword arguments to pass to make_adsorbate_structures()

        Returns
        -------
        Flow
            The Flow for this process.
        """
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

        if self.preset:
            self.bulk_to_slabs_job.preset = self.preset
            self.slab_to_adsorbates_job.preset = self.preset
        if self.swaps:
            self.bulk_to_slabs_job.swaps = self.swaps
            self.slab_to_adsorbates_job.swaps = self.swaps

        bulk_to_slabs_job = self.bulk_to_slabs_job.make(
            atoms, max_slabs=max_slabs, **slabgen_kwargs
        )
        jobs.append(bulk_to_slabs_job)

        if stable_slab:
            find_stable_slab_job = _get_slab_stability(
                bulk_static_job.output, bulk_to_slabs_job.output["all_outputs"]
            )
            slab_to_adsorbates_job = self.slab_to_adsorbates_job.make(
                find_stable_slab_job.output["stable_slab"]["atoms"],
                adsorbate,
                **make_ads_kwargs
            )
            jobs += [find_stable_slab_job, slab_to_adsorbates_job]
        else:
            slab_to_adsorbates_job = self.slab_to_adsorbates_job.make(
                bulk_to_slabs_job.output["all_atoms"], adsorbate, **make_ads_kwargs
            )
            jobs.append(slab_to_adsorbates_job)

        return Flow(jobs, output=slab_to_adsorbates_job.output, name=self.name)


@job
def _get_slab_stability(
    bulk_summary: Dict[str, Any], slab_summaries: Dict[str, Any]
) -> Dict[str, Any]:
    """
    A job that determine the most stable surface slab (based on cleavage energy) for
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
