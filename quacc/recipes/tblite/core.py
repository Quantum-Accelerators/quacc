"""
Core recipes for the tblite code
"""
from __future__ import annotations

from typing import Any, Dict

import covalent as ct
from ase.atoms import Atoms
from monty.dev import requires

from quacc.schemas.calc import summarize_opt_run, summarize_run
from quacc.util.calc import ideal_gas_thermo, run_ase_opt, run_ase_vib, run_calc

try:
    from tblite.ase import TBLite
except ImportError:
    TBLite = None


@ct.electron
@requires(
    TBLite,
    "tblite must be installed. Try pip install tblite[ase]",
)
def StaticJob(
    atoms: Atoms, method: str = "GFN2-xTB", tblite_kwargs: Dict[str, Any] | None = None
) -> Dict[str, Any]:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        .Atoms object
    method
        GFN1-xTB, GFN2-xTB, and IPEA1-xTB.
    tblite_kwargs
        Dictionary of custom kwargs for the tblite calculator.

    Returns
    -------
    Dict
        Summary of the run.
    """
    if tblite_kwargs is None:
        tblite_kwargs = {}

    atoms.calc = TBLite(method=method, **tblite_kwargs)
    new_atoms = run_calc(atoms)
    summary = summarize_run(new_atoms, input_atoms=atoms)

    return summary


@ct.electron
@requires(
    TBLite,
    "tblite must be installed. Try pip install tblite[ase]",
)
def RelaxJob(
    atoms: Atoms,
    method: str = "GFN2-xTB",
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "FIRE",
    tblite_kwargs: Dict[str, Any] | None = None,
    opt_kwargs: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """
    Function to relax a structure.

    Parameters
    ----------
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        .Optimizer class to use for the relaxation.
    tblite_kwargs
        Dictionary of custom kwargs for the tblite calculator.
    opt_kwargs
        Dictionary of kwargs for the optimizer.

    Returns
    -------
    Dict
        Summary of the run.
    """

    if tblite_kwargs is None:
        tblite_kwargs = {}
    if opt_kwargs is None:
        opt_kwargs = {}

    atoms.calc = TBLite(method=method, **tblite_kwargs)
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
@requires(TBLite, "tblite must be installed. Try pip install tblite[ase]")
def ThermoJob(
    atoms: Atoms,
    method: str = "GFN2-xTB",
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    xtb_kwargs: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """
    Function to run a frequency job and calculate thermochemistry.

    Parameters
    ----------
    atoms
        .Atoms object
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF.
    energy
        Potential energy in eV. If 0, then the output is just the correction.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    xtb_kwargs
        Dictionary of custom kwargs for the xTB calculator.

    Returns
    -------
    Dict
        Summary of the run.
    """
    if xtb_kwargs is None:
        xtb_kwargs = {}

    atoms.calc = TBLite(method=method, **xtb_kwargs)
    vibrations = run_ase_vib(atoms)
    thermo_summary = ideal_gas_thermo(
        atoms,
        vibrations.get_frequencies(),
        temperature=temperature,
        pressure=pressure,
        energy=energy,
    )

    return thermo_summary
