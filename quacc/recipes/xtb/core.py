"""
Core recipes for the xTB code
"""
from dataclasses import dataclass, field
from typing import Any, Dict

from ase.atoms import Atoms
from ase.optimize import FIRE
from ase.optimize.optimize import Optimizer
from jobflow import Maker, job
from monty.dev import requires

try:
    from xtb.ase.calculator import XTB
except (ModuleNotFoundError, ImportError):
    XTB = None
from quacc.schemas.calc import summarize_run


@dataclass
class StaticJob(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF.
    xtb_kwargs
        Dictionary of custom kwargs for the xTB calculator.
    """

    name: str = "xTB-Static"
    method: str = "GFN2-xTB"
    xtb_kwargs: Dict[str, Any] = field(default_factory=dict)

    @job
    @requires(
        XTB, "xTB-python must be installed. Try conda install -c conda-forge xtb-python"
    )
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
        atoms.calc = XTB(method=self.method, **self.xtb_kwargs)
        atoms.get_potential_energy()
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary


@dataclass
class RelaxJob(Maker):
    """
    Class to relax a structure.

    Parameters
    ----------
    name
        Name of the job.
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF.
    xtb_kwargs
        Dictionary of custom kwargs for the xTB calculator.
    optimizer
        .Optimizer class to use for the relaxation.
    fmax
        Tolerance for the force convergence (in eV/A).
    opt_kwargs
        Dictionary of kwargs for the optimizer.
    """

    name: str = "xTB-Relax"
    method: str = "GFN2-xTB"
    xtb_kwargs: Dict[str, Any] = field(default_factory=dict)
    optimizer: Optimizer = FIRE
    fmax: float = 0.03
    opt_kwargs: Dict[str, Any] = field(default_factory=dict)

    @job
    @requires(
        XTB, "xTB-python must be installed. Try conda install -c conda-forge xtb-python"
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
        # We always want to save the logfile and trajectory, so we will set some default
        # values if not specified by the user (and then remove them from the **opt_kwargs)
        logfile = self.opt_kwargs.get("logfile") or "opt.log"
        trajectory = self.opt_kwargs.get("trajectory") or "opt.traj"
        self.opt_kwargs.pop("logfile", None)
        self.opt_kwargs.pop("trajectory", None)

        atoms.calc = XTB(method=self.method, **self.xtb_kwargs)
        dyn = self.optimizer(
            atoms, logfile=logfile, trajectory=trajectory, **self.opt_kwargs
        )
        dyn.run(fmax=self.fmax)
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary
