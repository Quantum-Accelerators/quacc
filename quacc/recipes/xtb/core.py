"""
Core recipes for the xTB code
"""
from __future__ import annotations

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
from quacc.util.calc import calculate_thermo, run_ase_opt, run_ase_vib, run_calc


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
        summary = summarize_run(
            new_atoms, input_atoms=atoms, additional_fields={"name": self.name}
        )

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
    fmax
        Tolerance for the force convergence (in eV/A).
    optimizer
        .Optimizer class to use for the relaxation.
    xtb_kwargs
        Dictionary of custom kwargs for the xTB calculator.
    opt_kwargs
        Dictionary of kwargs for the optimizer.
    """

    name: str = "xTB-Relax"
    method: str = "GFN2-xTB"
    fmax: float = 0.01
    optimizer: Optimizer = FIRE
    xtb_kwargs: Dict[str, Any] = field(default_factory=dict)
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
        summary = summarize_run(
            new_atoms, input_atoms=atoms, additional_fields={"name": self.name}
        )

        return summary


@dataclass
class ThermoJob(Maker):
    """
    Class to calculate thermochemistry.

    Parameters
    ----------
    name
        Name of the job.
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    energy
        Potential energy in eV. If 0 eV, then the thermochemical correction is computed.
    geometry
        Monatomic, linear, or nonlinear. Will try to determine automatically if None.
    symmetry_number
        Rotational symmetry number.
    xtb_kwargs
        Dictionary of custom kwargs for the xTB calculator.
    """

    name: str = "xTB-Relax"
    method: str = "GFN2-xTB"
    temperature: float = 298.15
    pressure: float = 1.0
    energy: float = 0.0
    geometry: str = None
    symmetry_number: int = 1
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
            .Atoms object

        Returns
        -------
        Dict
            Summary of the thermochemistry.
        """
        atoms.calc = XTB(method=self.method, **self.xtb_kwargs)
        vibrations = run_ase_vib(atoms)
        thermo_summary = calculate_thermo(
            vibrations,
            temperature=self.temperature,
            pressure=self.pressure,
            energy=self.energy,
            geometry=self.geometry,
            symmetry_number=self.symmetry_number,
        )

        return thermo_summary
