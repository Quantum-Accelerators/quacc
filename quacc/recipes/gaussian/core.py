"""Core recipes for Gaussian"""
from __future__ import annotations

import multiprocessing
from dataclasses import dataclass, field
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.gaussian import Gaussian
from jobflow import Maker, job

from quacc.schemas.cclib import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc

LOG_FILE = Gaussian().label + ".log"
GEOM_FILE = LOG_FILE


@dataclass
class StaticJob(Maker):
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
    pop
        Type of population analysis to perform, if any
    molden
        Whether to write a molden file for orbital visualization
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "Gaussian-Static"
    xc: str = "wb97x-d"
    basis: str = "def2-tzvp"
    pop: str = "hirshfeld"
    molden: bool = True
    swaps: Dict[str, Any] = field(default_factory=dict)

    @job
    def make(
        self, atoms: Atoms, charge: int = None, mult: int = None
    ) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object
        charge
            Charge of the system. If None, this is determined from the sum of
            atoms.get_initial_charges().
        mult
            Multiplicity of the system. If None, this is determined from 1+ the sum
            of atoms.get_initial_magnetic_moments().

        Returns
        -------
        Dict
            Summary of the run.
        """
        defaults = {
            "mem": "16GB",
            "chk": "Gaussian.chk",
            "nprocshared": multiprocessing.cpu_count(),
            "xc": self.xc,
            "basis": self.basis,
            "charge": charge,
            "mult": mult,
            "sp": "",
            "scf": ["maxcycle=250", "xqc"],
            "integral": "ultrafine",
            "nosymmetry": "",
            "pop": self.pop,
            "gfinput": "" if self.molden else None,
            "ioplist": ["6/7=3"] if self.molden else None,
        }
        flags = merge_dicts(defaults, self.swaps, remove_none=True)

        atoms.calc = Gaussian(**flags)
        atoms = run_calc(atoms, geom_file=GEOM_FILE)
        summary = summarize_run(atoms, LOG_FILE, additional_fields={"name": self.name})

        return summary


@dataclass
class RelaxJob(Maker):
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
    freq
        If a requency calculation should be carried out.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "Gaussian-Relax"
    xc: str = "wb97x-d"
    basis: str = "def2-tzvp"
    freq: bool = False
    swaps: Dict[str, Any] = field(default_factory=dict)

    @job
    def make(
        self, atoms: Atoms, charge: int = None, mult: int = None
    ) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object
        charge
            Charge of the system. If None, this is determined from the sum of
            atoms.get_initial_charges().
        mult
            Multiplicity of the system. If None, this is determined from 1+ the sum
            of atoms.get_initial_magnetic_moments().

        Returns
        -------
        Dict
            Summary of the run.
        """
        defaults = {
            "mem": "16GB",
            "chk": "Gaussian.chk",
            "nprocshared": multiprocessing.cpu_count(),
            "xc": self.xc,
            "basis": self.basis,
            "charge": charge,
            "mult": mult,
            "opt": "",
            "scf": ["maxcycle=250", "xqc"],
            "integral": "ultrafine",
            "nosymmetry": "",
            "freq": "" if self.freq else None,
        }
        flags = merge_dicts(defaults, self.swaps, remove_none=True)

        atoms.calc = Gaussian(**flags)
        atoms = run_calc(atoms, geom_file=GEOM_FILE)
        summary = summarize_run(atoms, LOG_FILE, additional_fields={"name": self.name})

        return summary
