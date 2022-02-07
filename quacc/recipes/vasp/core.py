from dataclasses import dataclass
from typing import Any, Dict

from jobflow import Maker, job

from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run
from quacc.util.calc import run_calc
from quacc.util.json import unjsonify


@dataclass
class RelaxMaker(Maker):
    """
    Class to relax a structure.

    Parameters
    ----------
    name
        Name of the job.
    volume_relax
        True if a volume relaxation (ISIF = 3) should be performed.
        False if only the positions (ISIF = 2) should be updated.
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "Relax"
    volume_relax: bool = True
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

        if self.volume_relax:
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
class StaticMaker(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "Static"
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
            "ediff": 1e-6,
            "ismear": -5,
            "isym": 2,
            "laechg": True,
            "lcharg": True,
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
