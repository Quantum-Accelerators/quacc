"""
Core recipes for the tblite code
"""
from __future__ import annotations

import covalent as ct
from ase import Atoms
from monty.dev import requires

from quacc.schemas.ase import (
    summarize_opt_run,
    summarize_run,
    summarize_thermo_run,
    summarize_vib_run,
)
from quacc.util.atoms import copy_atoms
from quacc.util.calc import run_ase_opt, run_ase_vib, run_calc
from quacc.util.thermo import ideal_gas

try:
    from tblite.ase import TBLite
except ImportError:
    TBLite = None


@ct.electron
@requires(
    TBLite,
    "tblite must be installed. Try pip install tblite[ase]",
)
def static_job(
    atoms: Atoms, method: str = "GFN2-xTB", tblite_kwargs: dict | None = None
) -> dict:
    """
    Carry out a single-point calculation.

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
    summary
        Summary of the calculation.
    """
    tblite_kwargs = tblite_kwargs or {}
    input_atoms = copy_atoms(atoms)

    atoms.calc = TBLite(method=method, **tblite_kwargs)
    atoms = run_calc(atoms)
    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        additional_fields={"name": "TBLite Static"},
    )


@ct.electron
@requires(
    TBLite,
    "tblite must be installed. Try pip install tblite[ase]",
)
def relax_job(
    atoms: Atoms,
    method: str = "GFN2-xTB",
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "FIRE",
    tblite_kwargs: dict | None = None,
    opt_kwargs: dict | None = None,
) -> dict:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        .Atoms object
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
    summary
        Summary of the calculation.
    """

    tblite_kwargs = tblite_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    atoms.calc = TBLite(method=method, **tblite_kwargs)
    traj = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
    )

    return summarize_opt_run(
        traj, atoms.calc.parameters, additional_fields={"name": "TBLite Relax"}
    )


@requires(TBLite, "tblite must be installed. Try pip install tblite[ase]")
@ct.electron
def freq_job(
    atoms: Atoms,
    method: str = "GFN2-xTB",
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    xtb_kwargs: dict | None = None,
) -> dict:
    """
    Run a frequency job and calculate thermochemistry.

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
        dictionary of custom kwargs for the xTB calculator.

    Returns
    -------
    thermo_summary
        Summary of the thermochemistry.
    """

    xtb_kwargs = xtb_kwargs or {}
    input_atoms = copy_atoms(atoms)

    atoms.calc = TBLite(method=method, **xtb_kwargs)
    vibrations = run_ase_vib(atoms)

    igt = ideal_gas(input_atoms, vibrations.get_frequencies(), energy=energy)

    return {
        "vib": summarize_vib_run(
            vibrations, additional_fields={"name": "TBLite Vibrations"}
        ),
        "thermo": summarize_thermo_run(
            igt,
            temperature=temperature,
            pressure=pressure,
            additional_fields={"name": "TBLite Thermo"},
        ),
    }
