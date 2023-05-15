"""Core recipes for Gaussian"""
from __future__ import annotations

import multiprocessing
from typing import Any

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.gaussian import Gaussian

from quacc.schemas.cclib import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import merge_dicts

LOG_FILE = Gaussian().label + ".log"
GEOM_FILE = LOG_FILE


@ct.electron
def static_job(
    atoms: Atoms,
    charge: int | None = None,
    mult: int | None = None,
    xc: str = "wb97x-d",
    basis: str = "def2-tzvp",
    pop: str = "hirshfeld",
    write_molden: bool = True,
    swaps: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Function to carry out a single-point calculation.

    Parameters
    ----------
    atoms
        .Atoms object
    charge
        Charge of the system. If None, this is determined from the sum of
        atoms.get_initial_charges().
    mult
        Multiplicity of the system. If None, this is determined from 1+ the sum
        of atoms.get_initial_magnetic_moments().
    xc
        Exchange-correlation functional
    basis
        Basis set
    pop
        Type of population analysis to perform, if any
    write_molden
        Whether to write a molden file for orbital visualization
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    summary
        Dictionary of results from the calculation.
    """

    swaps = swaps or {}

    defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": multiprocessing.cpu_count(),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": mult,
        "sp": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "pop": pop,
        "gfinput": "" if write_molden else None,
        "ioplist": ["6/7=3"] if write_molden else None,
    }
    flags = merge_dicts(defaults, swaps, remove_none=True)

    atoms.calc = Gaussian(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE)
    summary = summarize_run(atoms, LOG_FILE)

    return summary


@ct.electron
def relax_job(
    atoms: Atoms,
    charge: int = None,
    mult: int = None,
    xc: str = "wb97x-d",
    basis: str = "def2-tzvp",
    freq: bool = False,
    swaps: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Function to carry out a geometry optimization.

    Parameters
    ----------
    atoms
        .Atoms object
    charge
        Charge of the system. If None, this is determined from the sum of
        atoms.get_initial_charges().
    mult
        Multiplicity of the system. If None, this is determined from 1+ the sum
        of atoms.get_initial_magnetic_moments().
    xc
        Exchange-correlation functional
    basis
        Basis set
    freq
        If a requency calculation should be carried out.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    summary
        Dictionary of results from the calculation.
    """

    swaps = swaps or {}

    defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": multiprocessing.cpu_count(),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": mult,
        "opt": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "freq": "" if freq else None,
    }
    flags = merge_dicts(defaults, swaps, remove_none=True)

    atoms.calc = Gaussian(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE)
    summary = summarize_run(atoms, LOG_FILE)

    return summary
