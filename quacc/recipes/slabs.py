from typing import Optional, Dict
from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run
from quacc.util.slabs import make_max_slabs_from_bulk, make_adsorbate_structures
from ase.io.jsonio import encode, decode
from jobflow import job, Flow, Response, Maker
from dataclasses import dataclass


@dataclass
class SlabRelaxMaker(Maker):
    name: str = "SlabRelax"
    preset: Optional[str] = None
    npar: int = 1
    kpar: int = 1

    @job
    def make(self, atoms_json: str, **kwargs) -> Dict:
        atoms = decode(atoms_json)
        flags = {
            "auto_dipole": True,
            "ediff": 1e-5,
            "ediffg": -0.02,
            "isif": 2,
            "ibrion": 2,
            "ismear": 0,
            "isym": 0,
            "kpar": self.kpar,
            "lcharg": False,
            "lwave": False,
            "ncore": self.ncore,
            "nsw": 200,
            "sigma": 0.05,
        }
        for k, v in kwargs.items():
            flags[k] = v

        atoms = SmartVasp(atoms, preset=self.preset, **kwargs)
        atoms.get_potential_energy()
        summary = summarize_run(atoms)

        return summary


@dataclass
class SlabStaticMaker(Maker):
    name: str = "SlabStatic"
    preset: Optional[str] = None
    npar: int = 1
    kpar: int = 1

    @job
    def make(self, atoms_json: str, **kwargs) -> Dict:
        atoms = decode(atoms_json)
        flags = {
            "auto_dipole": True,
            "ediff": 1e-6,
            "ismear": -5,
            "isym": 2,
            "kpar": self.kpar,
            "laechg": True,
            "lcharg": True,
            "lvhar": True,
            "lwave": True,
            "ncore": self.ncore,
            "nedos": 5001,
            "nsw": 0,
            "sigma": 0.05,
        }
        for k, v in kwargs.items():
            flags[k] = v

        atoms = SmartVasp(atoms, preset=self.preset, **kwargs)
        atoms.get_potential_energy()
        summary = summarize_run(atoms)

        return summary


@dataclass
class BulkToSlabMaker(Maker):
    name: str = "BulkToSlab"
    preset: Optional[str] = None
    npar: int = 1
    kpar: int = 1

    def make(
        self, atoms_json: str, max_slabs: Optional[int] = None, **slab_kwargs
    ) -> Response:
        atoms = decode(atoms_json)
        slabs = make_max_slabs_from_bulk(atoms, max_slabs=max_slabs, **slab_kwargs)
        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = SlabRelaxMaker(
                preset=self.preset, npar=self.npar, kpar=self.kpar
            ).make(encode(slab))
            jobs.append(relax_job)
            outputs.append(relax_job.output)

            static_job = SlabStaticMaker(
                preset=self.preset, npar=self.npar, kpar=self.kpar
            ).make(relax_job.output["atoms"])
            jobs.append(static_job)
            outputs.append(static_job.output)

        return Response(replace=Flow(jobs))


@dataclass
class SlabToAdsSlabMaker(Maker):
    name: str = "SlabToAdsSlab"
    preset: Optional[str] = None
    npar: int = 1
    kpar: int = 1

    def make(self, atoms_json: str, adsorbate_json: str, **slab_ads_kwargs) -> Response:
        atoms = decode(atoms_json)
        adsorbate = decode(adsorbate_json)

        slabs = make_adsorbate_structures(atoms, adsorbate, **slab_ads_kwargs)
        if slabs is None:
            return None

        jobs = []
        outputs = []
        for slab in slabs:
            relax_job = SlabRelaxMaker(
                preset=self.preset, npar=self.npar, kpar=self.kpar
            ).make(encode(slab))
            jobs.append(relax_job)
            outputs.append(relax_job.output)

            static_job = SlabStaticMaker(
                preset=self.preset, npar=self.npar, kpar=self.kpar
            ).make(relax_job.output["atoms"])
            jobs.append(static_job)
            outputs.append(static_job.output)

        return Response(replace=Flow(jobs))
