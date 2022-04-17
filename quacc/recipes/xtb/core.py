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
from quacc.util.calc import run_ase_opt, run_calc


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
        new_atoms = run_calc(atoms)
        summary = summarize_run(new_atoms, atoms, additional_fields={"name": self.name})

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
    fmax
        Tolerance for the force convergence (in eV/A).
    optimizer
        .Optimizer class to use for the relaxation.
    opt_kwargs
        Dictionary of kwargs for the optimizer.
    """

    name: str = "xTB-Relax"
    method: str = "GFN2-xTB"
    xtb_kwargs: Dict[str, Any] = field(default_factory=dict)
    fmax: float = 0.01
    optimizer: Optimizer = FIRE
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
        atoms.calc = XTB(method=self.method, **self.xtb_kwargs)
        new_atoms = run_ase_opt(
            atoms, fmax=self.fmax, optimizer=self.optimizer, opt_kwargs=self.opt_kwargs
        )
        summary = summarize_run(new_atoms, atoms, additional_fields={"name": self.name})

        return summary
