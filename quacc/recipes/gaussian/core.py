"""Core recipes for Gaussian"""
from __future__ import annotations

import multiprocessing

import covalent as ct
from ase import Atoms
from ase.calculators.gaussian import Gaussian

from quacc.schemas.cclib import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import merge_dicts

LOG_FILE = f"{Gaussian().label}.log"
GEOM_FILE = LOG_FILE


@ct.electron
def static_job(
    atoms: Atoms,
    charge: int = None,
    mult: int = None,
    xc: str = "wb97x-d",
    basis: str = "def2-tzvp",
    pop: str = "hirshfeld",
    write_molden: bool = True,
    swaps: dict = None,
) -> dict:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
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
    dict
        Dictionary of results from quacc.schemas.cclib.summarize_run
    """

    swaps = swaps or {}

    defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": multiprocessing.cpu_count(),
        "xc": xc,
        "basis": basis,
        "charge": charge or round(sum(atoms.get_initial_charges())),
        "mult": mult or round(1 + sum(atoms.get_initial_magnetic_moments())),
        "sp": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "pop": pop,
        "gfinput": "" if write_molden else None,
        "ioplist": ["6/7=3", "2/9=2000"]
        if write_molden
        else ["2/9=2000"],  # see ASE issue #660
    }
    flags = merge_dicts(defaults, swaps)

    atoms.calc = Gaussian(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE)

    return summarize_run(atoms, LOG_FILE, additional_fields={"name": "Gaussian Static"})


@ct.electron
def relax_job(
    atoms: Atoms,
    charge: int = None,
    mult: int = None,
    xc: str = "wb97x-d",
    basis: str = "def2-tzvp",
    freq: bool = False,
    swaps: dict = None,
) -> dict:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
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
        If a frequency calculation should be carried out.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.cclib.summarize_run
    """

    swaps = swaps or {}

    defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": multiprocessing.cpu_count(),
        "xc": xc,
        "basis": basis,
        "charge": charge or round(sum(atoms.get_initial_charges())),
        "mult": mult or round(1 + sum(atoms.get_initial_magnetic_moments())),
        "opt": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "freq": "" if freq else None,
        "ioplist": ["2/9=2000"],  # ASE issue #660
    }
    flags = merge_dicts(defaults, swaps)

    atoms.calc = Gaussian(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE)

    return summarize_run(atoms, LOG_FILE, additional_fields={"name": "Gaussian Relax"})
