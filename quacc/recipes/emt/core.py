from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.emt import EMT
from ase.optimize import FIRE
from ase.optimize.optimize import Optimizer
from jobflow import Maker, job

from quacc.schemas.calc import summarize_run

# NOTE: This set of recipes is mainly for demonstration purposes


@dataclass
class StaticMaker(Maker):
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
        atoms.calc = EMT(asap_cutoff=self.asap_cutoff)
        atoms.get_potential_energy()
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary


@dataclass
class RelaxMaker(Maker):
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
    opt_kwargs: Dict[str, Any] = None

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
        opt_kwargs = self.opt_kwargs or {}

        # We always want to save the logfile and trajectory, so we will set some default
        # values if not specified by the user (and then remove them from the **opt_kwargs)
        logfile = opt_kwargs.get("logfile", None) or "opt.log"
        trajectory = opt_kwargs.get("trajectory", None) or "opt.traj"
        opt_kwargs.pop("logfile", None)
        opt_kwargs.pop("trajectory", None)

        atoms.calc = EMT(asap_cutoff=self.asap_cutoff)
        dyn = self.optimizer(atoms,
                             logfile=logfile,
                             trajectory=trajectory,
                             **opt_kwargs)
        dyn.run(fmax=self.fmax)
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary
