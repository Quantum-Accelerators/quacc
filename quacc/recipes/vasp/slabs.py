from dataclasses import dataclass
from typing import Any, Dict

from jobflow import Flow, Maker, Response, job

from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run
from quacc.util.calc import run_calc
from quacc.util.json import jsonify, unjsonify
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

    name: str = "SlabRelax"
    preset: None | str = None
    swaps: Dict[str, Any] = None

    @job
    def make(self, atoms_json: str) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms_json
            Encoded .Atoms object

        Returns
        -------
        Dict:
            Summary of the run.
        """
        atoms = unjsonify(atoms_json)
        swaps = self.swaps or {}
        flags = {
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
        for k, v in swaps.items():
            flags[k] = v

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

    name: str = "SlabStatic"
    preset: None | str = None
    swaps: Dict[str, Any] = None

    @job
    def make(self, atoms_json: str) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms_json
            Encoded .Atoms object

        Returns
        -------
        Dict
            Summary of the run.
        """
        atoms = unjsonify(atoms_json)
        swaps = self.swaps or {}
        flags = {
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
        for k, v in swaps.items():
            flags[k] = v

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
    slab_relax_maker
        Maker to use for the SlabRelax job.
    slab_static_maker
        Default to use for the SlabStatic job.
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "BulkToSlab"
    slab_relax_maker: Maker = SlabRelaxMaker()
    slab_static_maker: Maker = SlabStaticMaker()
    preset: None | str = None
    swaps: Dict[str, Any] = None

    @job
    def make(
        self,
        atoms_json: str,
        max_slabs: None | int = None,
        slabgen_kwargs: Dict[str, Any] = None,
    ) -> Response:
        """
        Make the run.

        Parameters
        ----------
        atoms_json
            Encoded .Atoms object
        max_slabs
            Maximum number of slabs to make
        slabgen_kwargs
            Additional keyword arguments to pass to make_max_slabs_from_bulk()

        Returns
        -------
        Response
            A Flow of relaxation and static jobs for the generated slabs.
        """
        atoms = unjsonify(atoms_json)
        slabgen_kwargs = slabgen_kwargs or {}
        self.slab_relax_maker.preset = self.slab_relax_maker.preset or self.preset
        self.slab_static_maker.preset = self.slab_static_maker.preset or self.preset
        self.slab_relax_maker.swaps = self.slab_relax_maker.swaps or self.swaps
        self.slab_static_maker.swaps = self.slab_static_maker.swaps or self.swaps

        slabs = make_max_slabs_from_bulk(atoms, max_slabs=max_slabs, **slabgen_kwargs)
        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = self.slab_relax_maker.make(jsonify(slab))
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
    name:
        Name of the job.
    slab_relax_maker
        Maker to use for the SlabRelax job.
    slab_static_maker
        Default to use for the SlabStatic job.
    preset:
        Preset to use.
    """

    name: str = "SlabToAdsSlab"
    slab_relax_maker: Maker = SlabRelaxMaker()
    slab_static_maker: Maker = SlabStaticMaker()
    preset: None | str = None
    swaps: Dict[str, Any] = None

    @job
    def make(
        self,
        atoms_json: str,
        adsorbate_json: str,
        slabgen_ads_kwargs: Dict[str, Any] = None,
    ) -> Response | None:
        """
        Make the run.

        Parameters
        ----------
        atoms_json
            Encoded .Atoms object for the structure.
        adsorbate_json
            Encoded .Atoms object for the adsorbate.
        slabgen_ads_kwargs
            Additional keyword arguments to pass to make_adsorbate_structures()

        Returns
        -------
        Response
            A Flow of relaxation and static jobs for the generated slabs with adsorbates.
        """
        atoms = unjsonify(atoms_json)
        adsorbate = unjsonify(adsorbate_json)
        slabgen_ads_kwargs = slabgen_ads_kwargs or {}
        self.slab_relax_maker.preset = self.slab_relax_maker.preset or self.preset
        self.slab_static_maker.preset = self.slab_static_maker.preset or self.preset
        self.slab_relax_maker.swaps = self.slab_relax_maker.swaps or self.swaps
        self.slab_static_maker.swaps = self.slab_static_maker.swaps or self.swaps

        slabs = make_adsorbate_structures(atoms, adsorbate, **slabgen_ads_kwargs)
        if slabs is None:
            return None

        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = self.slab_relax_maker.make(jsonify(slab))
            jobs.append(relax_job)
            outputs.append(relax_job.output)

            static_job = self.slab_static_maker.make(relax_job.output["atoms"])
            jobs.append(static_job)
            outputs.append(static_job.output)

        return Response(replace=Flow(jobs))
