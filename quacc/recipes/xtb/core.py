"""
Core recipes for the xTB code
"""
from __future__ import annotations

import warnings
from typing import Any, Dict

import covalent as ct
from ase.atoms import Atoms
from monty.dev import requires

try:
    from xtb.ase.calculator import XTB
except ImportError:
    XTB = None
from quacc.schemas.calc import summarize_opt_run, summarize_run
from quacc.util.calc import ideal_gas_thermo, run_ase_opt, run_ase_vib, run_calc


@ct.electron
@requires(XTB, "xTB-python must be installed. Try pip install xtb")
def StaticJob(
    atoms: Atoms, method: str = "GFN2-xTB", xtb_kwargs: Dict[str, Any] | None = None
) -> Dict[str, Any]:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        .Atoms object
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF.
    xtb_kwargs
        Dictionary of custom kwargs for the xTB calculator.

    Returns
    -------
    summary
        Summary of the calculation.
    """

    warnings.warn(
        "xTB-python will be deprecated in a future version of QuAcc. If possible, you should try using tblite instead.",
        DeprecationWarning,
    )

    xtb_kwargs = xtb_kwargs or {}

    atoms.calc = XTB(method=method, **xtb_kwargs)
    new_atoms = run_calc(atoms)
    summary = summarize_run(new_atoms, input_atoms=atoms)

    return summary


@ct.electron
@requires(XTB, "xTB-python must be installed. Try pip install xtb")
def RelaxJob(
    atoms: Atoms,
    method: str = "GFN2-xTB",
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "FIRE",
    xtb_kwargs: Dict[str, Any] | None = None,
    opt_kwargs: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """
    Function to relax a structure.

    Parameters
    ----------
    atoms
        .Atoms object
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        Name of ASE optimizer class to use for the relaxation.
    xtb_kwargs
        Dictionary of custom kwargs for the xTB calculator.
    opt_kwargs
        Dictionary of kwargs for the optimizer.

    Returns
    -------
    summary
        Summary of the calculation.
    """

    xtb_kwargs = xtb_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    atoms.calc = XTB(method=method, **xtb_kwargs)
    traj = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
    )
    summary = summarize_opt_run(traj, atoms.calc.parameters)

    return summary


@ct.electron
@requires(XTB, "xTB-python must be installed. Try pip install xtb")
def ThermoJob(
    atoms: Atoms,
    method: str = "GFN2-xTB",
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    xtb_kwargs: Dict[str, Any] | None = None,
):
    """
    Function to run a frequency job and calculate thermochemistry.

    Parameters
    ----------
    atoms
        .Atoms object
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    xtb_kwargs
        Dictionary of custom kwargs for the xTB calculator.

    Returns
    -------
    thermo_summary
        Summary of the thermochemistry.
    """

    xtb_kwargs = xtb_kwargs or {}

    atoms.calc = XTB(method=method, **xtb_kwargs)
    vibrations = run_ase_vib(atoms)
    thermo_summary = ideal_gas_thermo(
        atoms,
        vibrations.get_frequencies(),
        temperature=temperature,
        pressure=pressure,
        energy=energy,
    )

    return thermo_summary
