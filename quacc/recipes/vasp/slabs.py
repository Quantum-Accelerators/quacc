from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from jobflow import Flow, Maker, Response, job

from quacc.calculators.vasp import SmartVasp
from quacc.recipes.vasp.core import RelaxJob, StaticJob
from quacc.schemas.vasp import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc
from quacc.util.slabs import make_adsorbate_structures, make_max_slabs_from_bulk


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

        return Response(replace=Flow(jobs))  # , output=outputs))


@dataclass
class SlabToAdsSlabJob(Maker):
    """
    Class to convert a slab structure to one with adsorbates present,
    along with the relaxations and statics for the slab-adsorbate systems.

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

    name: str = "VASP-SlabToAdsSlab"
    preset: str = None
    swaps: Dict[str, Any] = None
    slab_relax_job: Maker | None = SlabRelaxJob()
    slab_static_job: Maker | None = SlabStaticJob()

    @job
    def make(self, atoms: Atoms, adsorbate: Atoms, **slabgen_ads_kwargs) -> Response:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object for the structure.
        adsorbate
            .Atoms object for the adsorbate.
        slabgen_ads_kwargs
            Additional keyword arguments to pass to make_adsorbate_structures()

        Returns
        -------
        Response
            A Flow of relaxation and static jobs for the generated slabs with adsorbates.
        """
        slabgen_ads_kwargs = slabgen_ads_kwargs or {}
        if self.preset:
            self.slab_static_job.preset = self.preset
            self.slab_relax_job.preset = self.preset
        if self.swaps:
            self.slab_static_job.swaps = self.swaps
            self.slab_relax_job.swaps = self.swaps

        slabs = make_adsorbate_structures(atoms, adsorbate, **slabgen_ads_kwargs)

        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = self.slab_relax_job.make(slab)
            jobs.append(relax_job)

            static_job = self.slab_static_job.make(relax_job.output["atoms"])
            jobs.append(static_job)
            outputs.append(static_job.output)

        return Response(replace=Flow(jobs))  # , output=outputs))


@dataclass
class BulktoAdsEnergyFlow(Maker):
    """
    Maker to get adsorption energies from a bulk structure.
    """

    name: str = "Bulk"
    preset: str = None
    bulk_relax_job: Maker | None = RelaxJob()
    bulk_static_job: Maker | None = StaticJob()
    bulk_to_slabs_job: Maker = BulkToSlabsJob()
    slab_to_ads_slab_job: Maker = SlabToAdsSlabJob()
    swaps: Dict[str, Any] = None

    def make(self, atoms: Atoms, adsorbate: Atoms):
        jobs = []

        if self.bulk_relax_job:
            if self.preset:
                bulk_relax_job.preset = self.preset
            if self.swaps:
                bulk_relax_job.swaps = self.swaps
            bulk_relax_job = self.bulk_relax_job.make(atoms)
            atoms = bulk_relax_job.output["atoms"]
            jobs.append(bulk_relax_job)

        if self.bulk_static_job:
            if self.preset:
                bulk_static_job.preset = self.preset
            if self.swaps:
                bulk_static_job.swaps = self.swaps
            bulk_static_job = self.bulk_static_job.make(atoms)
            atoms = bulk_static_job.output["atoms"]
            jobs.append(bulk_static_job)

        bulk_to_slabs_job = self.bulk_to_slabs_job.make(atoms)

        # bulk_to_slabs_job.output

        jobs += [bulk_to_slabs_job]
        return Flow(jobs)
