from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from jobflow import Flow, Maker, Response, job

from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc
from quacc.util.slabs import make_adsorbate_structures, make_max_slabs_from_bulk


@dataclass
class SlabRelaxMaker(Maker):
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
    preset: None | str = None
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
            "ediff": 1e-5,
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
class SlabStaticMaker(Maker):
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
    preset: None | str = None
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
            "ediff": 1e-6,
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
class BulkToSlabMaker(Maker):
    """
    Class to convert a bulk structure to a slab,
    along with the relaxations and statics for the slabs.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use. Applies to all jobs in the flow.
    slab_relax_maker
        Maker to use for the SlabRelax job.
    slab_static_maker
        Default to use for the SlabStatic job.
    swaps
        Dictionary of custom kwargs for the calculator.
        Applies to all jobs in the flow.
    """

    name: str = "VASP-BulkToSlab"
    preset: None | str = None
    slab_relax_maker: Maker = SlabRelaxMaker()
    slab_static_maker: Maker = SlabStaticMaker()
    swaps: Dict[str, Any] = None

    @job
    def make(self,
             atoms: Atoms,
             max_slabs: None | int = None,
             **slabgen_kwargs) -> Response:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object
        max_slabs
            Maximum number of slabs to make
        slabgen_kwargs
            Additional keyword arguments to pass to make_max_slabs_from_bulk()

        Returns
        -------
        Response
            A Flow of relaxation and static jobs for the generated slabs.
        """
        slabgen_kwargs = slabgen_kwargs or {}
        self.slab_relax_maker.preset = self.slab_relax_maker.preset or self.preset
        self.slab_static_maker.preset = self.slab_static_maker.preset or self.preset
        self.slab_relax_maker.swaps = self.slab_relax_maker.swaps or self.swaps
        self.slab_static_maker.swaps = self.slab_static_maker.swaps or self.swaps

        slabs = make_max_slabs_from_bulk(atoms,
                                         max_slabs=max_slabs,
                                         **slabgen_kwargs)
        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = self.slab_relax_maker.make(slab)
            jobs.append(relax_job)
            outputs.append(relax_job.output)

            static_job = self.slab_static_maker.make(relax_job.output["atoms"])
            jobs.append(static_job)
            outputs.append(static_job.output)

        return Response(replace=Flow(jobs))


@dataclass
class SlabToAdsSlabMaker(Maker):
    """
    Class to convert a slab structure to one with adsorbates present,
    along with the relaxations and statics for the slab-adsorbate systems.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use. Applies to all jobs in the flow.
    slab_relax_maker
        Maker to use for the SlabRelax job.
    slab_static_maker
        Default to use for the SlabStatic job.
    swaps
        Dictionary of custom kwargs for the calculator.
        Applies to all jobs in the flow.
    """

    name: str = "VASP-SlabToAdsSlab"
    preset: None | str = None
    swaps: Dict[str, Any] = None
    slab_relax_maker: Maker = SlabRelaxMaker()
    slab_static_maker: Maker = SlabStaticMaker()

    @job
    def make(self, atoms: Atoms, adsorbate: Atoms | str,
             **slabgen_ads_kwargs) -> Response:
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
        self.slab_relax_maker.preset = self.slab_relax_maker.preset or self.preset
        self.slab_static_maker.preset = self.slab_static_maker.preset or self.preset
        self.slab_relax_maker.swaps = self.slab_relax_maker.swaps or self.swaps
        self.slab_static_maker.swaps = self.slab_static_maker.swaps or self.swaps

        slabs = make_adsorbate_structures(atoms, adsorbate,
                                          **slabgen_ads_kwargs)

        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = self.slab_relax_maker.make(slab)
            jobs.append(relax_job)
            outputs.append(relax_job.output)

            static_job = self.slab_static_maker.make(relax_job.output["atoms"])
            jobs.append(static_job)
            outputs.append(static_job.output)

        return Response(replace=Flow(jobs))
