"""Core recipes for the tblite code"""
from __future__ import annotations

from typing import Literal

import covalent as ct
from ase import Atoms
from ase.optimize import FIRE
from monty.dev import requires

from quacc.schemas.ase import (
    DynSchema,
    RunSchema,
    ThermoSchema,
    VibSchema,
    summarize_dyn_run,
    summarize_run,
    summarize_thermo_run,
    summarize_vib_run,
)
from quacc.util.calc import run_ase_opt, run_ase_vib, run_calc
from quacc.util.thermo import ideal_gas

try:
    from tblite.ase import TBLite
except ImportError:
    TBLite = None


@ct.electron
@requires(TBLite, "tblite must be installed. Try pip install tblite[ase]")
def static_job(
    atoms: Atoms | dict,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    method
        GFN1-xTB, GFN2-xTB, and IPEA1-xTB.
    calc_swaps
        Dictionary of custom kwargs for the tblite calculator.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from quacc.schemas.ase.summarize_run
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    calc_swaps = calc_swaps or {}

    atoms.calc = TBLite(method=method, **calc_swaps)
    final_atoms = run_calc(atoms, copy_files=copy_files)
    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "TBLite Static"},
    )


@ct.electron
@requires(TBLite, "tblite must be installed. Try pip install tblite[ase]")
def relax_job(
    atoms: Atoms | dict,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> DynSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB.
    calc_swaps
        Dictionary of custom kwargs for the tblite calculator.
    opt_swaps
        Dictionary of custom kwargs for run_ase_opt
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    DynSchema
        Dictionary of results from quacc.schemas.ase.summarize_dyn_run
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    calc_swaps = calc_swaps or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = opt_defaults | opt_swaps

    atoms.calc = TBLite(method=method, **calc_swaps)
    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_flags)

    return summarize_dyn_run(dyn, additional_fields={"name": "TBLite Relax"})


@ct.electron
@requires(TBLite, "tblite must be installed. Try pip install tblite[ase]")
def freq_job(
    atoms: Atoms | dict,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    calc_swaps: dict | None = None,
    vib_kwargs: dict | None = None,
    copy_files: list[str] | None = None,
) -> dict[Literal["vib", "thermo"], VibSchema | ThermoSchema]:
    """
    Run a frequency job and calculate thermochemistry.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF.
    energy
        Potential energy in eV. If 0, then the output is just the correction.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    calc_swaps
        dictionary of custom kwargs for the xTB calculator.
    vib_kwargs
        dictionary of custom kwargs for the Vibrations object.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.ase.summarize_vib_run and
        quacc.schemas.ase.summarize_thermo_run
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    calc_swaps = calc_swaps or {}
    vib_kwargs = vib_kwargs or {}

    atoms.calc = TBLite(method=method, **calc_swaps)
    vibrations = run_ase_vib(atoms, vib_kwargs=vib_kwargs, copy_files=copy_files)

    igt = ideal_gas(atoms, vibrations.get_frequencies(), energy=energy)

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
