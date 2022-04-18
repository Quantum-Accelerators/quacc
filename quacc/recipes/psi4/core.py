"""Core recipes for Psi4"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.psi4 import Psi4
from jobflow import Maker, job
from monty.dev import requires

try:
    import psi4
except:
    psi4 = None
from quacc.schemas.calc import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc


@dataclass
class StaticMaker(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    method
        The level of theory to use.
    basis
        Basis set
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "psi4-Static"
    method: str = "wb97x-v"
    basis: str = "def2-tzvp"
    swaps: Dict[str, Any] = None

    @job
    @requires(psi4, "Psi4 be installed. Try conda install -c psi4 psi4")
    def make(
        self, atoms: Atoms, charge: int = None, mult: int = None
    ) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object`
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
        swaps = self.swaps or {}
        defaults = {
            "mem": "16GB",
            "num_thread": "max",
            "method": self.method,
            "basis": "def2-tzvp",
            "charge": charge if charge else round(sum(atoms.get_initial_charges())),
            "multiplicity": mult
            if mult
            else round(1 + sum(atoms.get_initial_magnetic_moments())),
        }
        flags = merge_dicts(defaults, swaps, remove_none=True)

        atoms.calc = Psi4(**flags)
        atoms = run_calc(atoms)
        summary = summarize_run(atoms, ".log", additional_fields={"name": self.name})

        return summary
