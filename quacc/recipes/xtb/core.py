from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from ase.optimize import FIRE
from ase.optimize.optimize import Optimizer
from jobflow import Maker, job
from monty.dev import requires

try:
    from xtb.ase.calculator import XTB
except:
    XTB = None
from quacc.schemas.calc import summarize_run


@dataclass
class StaticMaker(Maker):
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
    xtb_kwargs: Dict[str, Any] = None

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
        xtb_kwargs = self.xtb_kwargs or {}

        atoms.calc = XTB(method=self.method, **xtb_kwargs)
        atoms.get_potential_energy()
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary


@dataclass
class RelaxMaker(Maker):
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
    optkwargs
        Dictionary of kwargs for the optimizer.
    """

    name: str = "xTB-Relax"
    method: str = "GFN2-xTB"
    xtb_kwargs: Dict[str, Any] = None
    optimizer: Optimizer = FIRE
    fmax: float = 0.03
    optkwargs: Dict[str, Any] = None

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
        xtb_kwargs = self.xtb_kwargs or {}
        optkwargs = self.optkwargs or {}
        logfile = optkwargs.get("logfile", None) or "opt.log"
        trajectory = optkwargs.get("trajectory", None) or "opt.traj"
        restart = optkwargs.get("restart", None) or "opt.pckl"

        atoms.calc = XTB(method=self.method, **xtb_kwargs)
        dyn = self.optimizer(
            atoms, logfile=logfile, trajectory=trajectory, restart=restart, **optkwargs
        )
        dyn.run(fmax=self.fmax)
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary
