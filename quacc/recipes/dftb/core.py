"""Core recipes for DFTB+"""
from dataclasses import dataclass, field
from shutil import which
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.dftb import Dftb
from jobflow import Maker, job
from monty.dev import requires

from quacc.schemas.calc import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc

DFTBPLUS_EXISTS = bool(which("dftb+"))


@dataclass
class StaticJob(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    Method
        Method to use. Accepts 'DFTB', 'GFN1-xTB', and 'GFN2-xTB'.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "DFTB-Static"
    method: str = "GFN2-xTB"
    swaps: Dict[str, Any] = field(default_factory=dict)

    @job
    @requires(
        DFTBPLUS_EXISTS,
        "DFTB+ must be installed. Try conda install -c conda-forge dftbplus",
    )
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
        defaults = {}
        if "xtb" in self.method.lower():
            defaults["Hamiltonian_"] = "xTB"
        if "gfn2-xtb" in self.method.lower():
            defaults["Hamiltonian_Method"] = "GFN2-xTB"
        elif "gfn1-xtb" in self.method.lower():
            defaults["Hamiltonian_Method"] = "GFN1-xTB"

        flags = merge_dicts(defaults, self.swaps, remove_none=True, auto_lowercase=False)

        atoms.calc = Dftb(**flags)
        atoms = run_calc(atoms)
        summary = summarize_run(
            atoms, "dftb.out", additional_fields={"name": self.name}
        )

        return summary
