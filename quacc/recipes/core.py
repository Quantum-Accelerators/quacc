from typing import Optional, Dict
from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run
from ase.io.jsonio import decode
from jobflow import job, Maker
from dataclasses import dataclass


@dataclass
class RelaxMaker(Maker):
    name: str = "Relax"
    preset: Optional[str] = None
    npar: int = 1
    kpar: int = 1

    @job
    def make(
        self, atoms_json: str, volume_relax: bool = True, **kwargs
    ) -> Dict:
        atoms = decode(atoms_json)
        if volume_relax:
            isif = 3
        else:
            isif = 2
        flags = {
            "ediff": 1e-5,
            "ediffg": -0.02,
            "isif": isif,
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
class StaticMaker(Maker):
    name: str = "Static"
    preset: Optional[str] = None
    npar: int = 1
    kpar: int = 1

    @job
    def make(self, atoms_json: str, slab: bool = False, **kwargs) -> Dict:
        atoms = decode(atoms_json)
        flags = {
            "ediff": 1e-6,
            "ismear": -5,
            "isym": 2,
            "kpar": self.kpar,
            "laechg": True,
            "lcharg": True,
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
