from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.gaussian import Gaussian
from jobflow import Maker, job

from quacc.schemas.cclib import summarize_run
from quacc.util.calc import run_calc


@dataclass
class StaticMaker(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    xc
        Exchange-correlation functional
    basis
        Basis set
    swaps
        Dictionary of custom kwargs for the calculator.
    pop
        Type of population analysis to perform, if any
    molden
        Whether to write a molden file for orbital visualization
    """

    name: str = "Gaussian-Static"
    xc: str = "wB97X-D"
    basis: str = "def2-TZVP"
    swaps: Dict[str, Any] = None
    pop: str = "hirshfeld"
    molden: bool = True

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
            "mem": "16GB",
            "xc": self.xc,
            "basis": self.basis,
            "sp": "",
            "scf": ["maxcycle=250", "xqc"],
            "integral": "ultrafine",
            "nosymmetry": "",
            "pop": self.pop,
            "gfinput": "" if self.molden else None,
            "ioplist": ["6/7=3"] if self.molden else None,
        }
        for k, v in swaps.items():
            flags[k] = v
        atoms.calc = Gaussian(**flags)
        atoms = run_calc(atoms)
        summary = summarize_run(atoms, ".log", additional_fields={"name": self.name})

        return summary


@dataclass
class RelaxMaker(Maker):
    """
    Class to carry out a geometry optimization.

    Parameters
    ----------
    name
        Name of the job.
    xc
        Exchange-correlation functional
    basis
        Basis set
    swaps
        Dictionary of custom kwargs for the calculator.
    freq
        If a requency calculation should be carried out.
    """

    name: str = "Gaussian-Relax"
    xc: str = "wB97X-D"
    basis: str = "def2-TZVP"
    swaps: Dict[str, Any] = None
    freq: bool = False

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
            "mem": "16GB",
            "xc": self.xc,
            "basis": self.basis,
            "opt": "",
            "scf": ["maxcycle=250", "xqc"],
            "integral": "ultrafine",
            "nosymmetry": "",
            "freq": "" if self.freq else None,
        }
        for k, v in swaps.items():
            flags[k] = v
        atoms.calc = Gaussian(**flags)
        atoms = run_calc(atoms)
        summary = summarize_run(atoms, ".log", additional_fields={"name": self.name})

        return summary
