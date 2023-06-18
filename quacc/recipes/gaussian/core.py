"""Core recipes for Gaussian"""
from __future__ import annotations

import multiprocessing

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.gaussian import Gaussian

from quacc.schemas.cclib import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import remove_dict_empties

LOG_FILE = f"{Gaussian().label}.log"
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
    swaps: dict | None = None,
) -> dict:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system. If None, this is determined from `atoms.charge`
    mult
        Multiplicity of the system. If None, this is determined from `atoms.spin_multiplicity`
    xc
        Exchange-correlation functional
    basis
        Basis set
    pop
        Type of population analysis to perform from `quacc.schemas.cclib.summarize_run`
    write_molden
        Whether to write a molden file for orbital visualization
    swaps
        Dictionary of custom kwargs for the calculator.
            defaults = {
                "mem": "16GB",
                "chk": "Gaussian.chk",
                "nprocshared": multiprocessing.cpu_count(),
                "xc": xc,
                "basis": basis,
                "charge": charge or getattr(atoms, "charge"),
                "mult": mult or getattr(atoms, "spin_multiplicity"),
                "sp": "",
                "scf": ["maxcycle=250", "xqc"],
                "integral": "ultrafine",
                "nosymmetry": "",
                "pop": pop,
                "gfinput": "" if write_molden else None,
                "ioplist": ["6/7=3"] if write_molden else [],
            }

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.cclib.summarize_run`
    """

    swaps = swaps or {}

    atoms.charge = charge or getattr(atoms, "charge")
    atoms.spin_multiplicity = mult or getattr(atoms, "spin_multiplicity")

    defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": multiprocessing.cpu_count(),
        "xc": xc,
        "basis": basis,
        "charge": atoms.charge,
        "mult": atoms.spin_multiplicity,
        "sp": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "pop": pop,
        "gfinput": "" if write_molden else None,
        "ioplist": ["6/7=3"] if write_molden else [],
    }
    flags = remove_dict_empties(defaults | swaps)

    atoms.calc = Gaussian(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE)

    return summarize_run(atoms, LOG_FILE, additional_fields={"name": "Gaussian Static"})


@ct.electron
def relax_job(
    atoms: Atoms,
    charge: int | None = None,
    mult: int | None = None,
    xc: str = "wb97x-d",
    basis: str = "def2-tzvp",
    freq: bool = False,
    swaps: dict | None = None,
) -> dict:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system. If None, this is determined from atoms.charge
    mult
        Multiplicity of the system. If None, this is determined from atoms.spin_multiplicity
    xc
        Exchange-correlation functional
    basis
        Basis set
    freq
        If a frequency calculation should be carried out.
    swaps
        Dictionary of custom kwargs for the calculator.
            defaults = {
                "mem": "16GB",
                "chk": "Gaussian.chk",
                "nprocshared": multiprocessing.cpu_count(),
                "xc": xc,
                "basis": basis,
                "charge": charge or getattr(atoms, "charge"),
                "mult": mult or getattr(atoms, "spin_multiplicity"),
                "opt": "",
                "scf": ["maxcycle=250", "xqc"],
                "integral": "ultrafine",
                "nosymmetry": "",
                "freq": "" if freq else None,
            }

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.cclib.summarize_run`
    """

    swaps = swaps or {}

    atoms.charge = charge or getattr(atoms, "charge")
    atoms.spin_multiplicity = mult or getattr(atoms, "spin_multiplicity")

    defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": multiprocessing.cpu_count(),
        "xc": xc,
        "basis": basis,
        "charge": atoms.charge,
        "mult": atoms.spin_multiplicity,
        "opt": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "freq": "" if freq else None,
    }
    flags = remove_dict_empties(defaults | swaps)

    atoms.calc = Gaussian(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE)

    return summarize_run(atoms, LOG_FILE, additional_fields={"name": "Gaussian Relax"})
