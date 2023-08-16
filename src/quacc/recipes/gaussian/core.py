"""Core recipes for Gaussian"""
from __future__ import annotations

import multiprocessing
from typing import TYPE_CHECKING

import covalent as ct
from ase.calculators.gaussian import Gaussian

from quacc.schemas.atoms import fetch_atoms
from quacc.schemas.cclib import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import get_parameters

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.cclib import cclibSchema

LOG_FILE = f"{Gaussian().label}.log"
GEOM_FILE = LOG_FILE


@ct.electron
def static_job(
    atoms: Atoms | dict,
    charge: int | None = None,
    multiplicity: int | None = None,
    xc: str = "wb97x-d",
    basis: str = "def2-tzvp",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> cclibSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    charge
        Charge of the system. If None, this is determined from the sum of
        `atoms.get_initial_charges().`
    multiplicity
        Multiplicity of the system. If None, this is determined from 1+ the sum
        of `atoms.get_initial_magnetic_moments()`.
    xc
        Exchange-correlation functional
    basis
        Basis set
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from `quacc.schemas.cclib.summarize_run`
    """
    atoms = fetch_atoms(atoms)

    charge = int(atoms.get_initial_charges().sum()) if charge is None else charge
    multiplicity = (
        int(1 + atoms.get_initial_magnetic_moments().sum())
        if multiplicity is None
        else multiplicity
    )

    defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": multiprocessing.cpu_count(),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": multiplicity,
        "sp": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "pop": "CM5",
        "gfinput": "",
        "ioplist": ["6/7=3", "2/9=2000"],  # see ASE issue #660
    }
    flags = get_parameters(defaults, calc_swaps)

    atoms.calc = Gaussian(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    return summarize_run(
        atoms,
        LOG_FILE,
        additional_fields={"name": "Gaussian Static"},
    )


@ct.electron
def relax_job(
    atoms: Atoms,
    charge: int | None = None,
    multiplicity: int | None = None,
    xc: str = "wb97x-d",
    basis: str = "def2-tzvp",
    freq: bool = False,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> cclibSchema:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    charge
        Charge of the system. If None, this is determined from the sum of
        `atoms.get_initial_charges()`.
    multiplicity
        Multiplicity of the system. If None, this is determined from 1+ the sum
        of `atoms.get_initial_magnetic_moments()`.
    xc
        Exchange-correlation functional
    basis
        Basis set
    freq
        If a frequency calculation should be carried out.
    calc_swaps
        Dictionary of custom kwargs for the calculator.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from `quacc.schemas.cclib.summarize_run`
    """
    atoms = fetch_atoms(atoms)

    charge = int(atoms.get_initial_charges().sum()) if charge is None else charge
    multiplicity = (
        int(1 + atoms.get_initial_magnetic_moments().sum())
        if multiplicity is None
        else multiplicity
    )

    defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": multiprocessing.cpu_count(),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": multiplicity,
        "opt": "",
        "pop": "CM5",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "freq": "" if freq else None,
        "ioplist": ["2/9=2000"],  # ASE issue #660
    }
    flags = get_parameters(defaults, calc_swaps)

    atoms.calc = Gaussian(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    return summarize_run(
        atoms,
        LOG_FILE,
        additional_fields={"name": "Gaussian Relax"},
    )
