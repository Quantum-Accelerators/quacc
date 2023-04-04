"""Core recipes for EMT"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.emt import EMT
from ase.io import read
from ase.optimize import FIRE
from jobflow import Maker, job

from quacc.schemas.calc import summarize_opt_run, summarize_run

# NOTE: This set of minimal recipes is mainly for demonstration purposes


@dataclass
class StaticJob(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    asap_cutoff
        If an ASAP-style cutoff should be used.
    """

    name: str = "EMT-Static"
    asap_cutoff: bool = False

    @job
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object`

        Returns
        -------
        Dict
            Summary of the run.
        """
        input_atoms = deepcopy(atoms)
        atoms.calc = EMT(asap_cutoff=self.asap_cutoff)
        atoms.get_potential_energy()
        summary = summarize_run(
            atoms, input_atoms=input_atoms, additional_fields={"name": self.name}
        )

        return summary


@dataclass
class RelaxJob(Maker):
    """
    Class to carry out a geometry optimization.

    Parameters
    ----------
    name
        Name of the job.
    asap_cutoff
        If an ASAP-style cutoff should be used.
    fmax
        Tolerance for the force convergence (in eV/A).
    opt_kwargs
        Dictionary of kwargs for the optimizer.
    """

    name: str = "EMT-Relax"
    asap_cutoff: bool = False
    fmax: float = 0.03
    opt_kwargs: Dict[str, Any] = field(default_factory=dict)

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
        atoms.calc = EMT(asap_cutoff=self.asap_cutoff)
        dyn = FIRE(atoms, logfile="opt.log", trajectory="opt.traj", **self.opt_kwargs)
        dyn.run(fmax=self.fmax)
        traj = read("opt.traj", index=":")
        summary = summarize_opt_run(
            traj, atoms.calc.parameters, additional_fields={"name": self.name}
        )

        return summary
