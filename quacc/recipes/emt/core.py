"""Core recipes for EMT"""
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.emt import EMT
from ase.optimize import FIRE
from ase.optimize.optimize import Optimizer
from jobflow import Maker, job

from quacc.schemas.calc import summarize_run
from quacc.util.calc import run_ase_opt, run_calc

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
            atoms, input_atoms, additional_fields={"name": self.name}
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
    optimizer
        .Optimizer class to use for the relaxation.
    fmax
        Tolerance for the force convergence (in eV/A).
    opt_kwargs
        Dictionary of kwargs for the optimizer.
    """

    name: str = "EMT-Relax"
    asap_cutoff: bool = False
    optimizer: Optimizer = FIRE
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
        input_atoms = deepcopy(atoms)
        atoms.calc = EMT(asap_cutoff=self.asap_cutoff)
        dyn = self.optimizer(
            atoms, logfile="opt.log", trajectory="opt.traj", **self.opt_kwargs
        )
        dyn.run(fmax=self.fmax)
        summary = summarize_run(
            atoms, input_atoms, additional_fields={"name": self.name}
        )

        return summary
