from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from jobflow import Maker, job

from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run
from quacc.util.calc import run_calc


@dataclass
class RelaxMaker(Maker):
    """
    Class to relax a structure.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    volume_relax
        True if a volume relaxation (ISIF = 3) should be performed.
        False if only the positions (ISIF = 2) should be updated.
    """

    name: str = "VASP-Relax"
    preset: None | str = None
    swaps: Dict[str, Any] = None
    volume_relax: bool = True

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
        flags = {
            "ediff": 1e-5,
            "ediffg": -0.02,
            "isif": 3 if self.volume_relax else 2,
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

    name: str = "VASP-Static"
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
