"""Core recipes for Gaussian"""
from __future__ import annotations

import multiprocessing
from typing import TYPE_CHECKING

from ase.calculators.gaussian import Gaussian

from quacc import job
from quacc.schemas import fetch_atoms
from quacc.schemas.cclib import summarize_run
from quacc.utils.atoms import get_charge, get_multiplicity
from quacc.utils.calc import run_calc
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.cclib import cclibSchema

LOG_FILE = f"{Gaussian().label}.log"
GEOM_FILE = LOG_FILE


@job
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
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
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
        Dictionary of custom kwargs for the calculator. Overrides the following
        defaults:

        ```python
        {
            "mem": "16GB",
            "chk": "Gaussian.chk",
            "nprocshared": multiprocessing.cpu_count(),
            "xc": xc,
            "basis": basis,
            "charge": get_charge(atoms) if charge is None else charge,
            "mult": get_multiplicity(atoms) if multiplicity is None else multiplicity,
            "sp": "",
            "scf": ["maxcycle=250", "xqc"],
            "integral": "ultrafine",
            "nosymmetry": "",
            "pop": "CM5",
            "gfinput": "",
            "ioplist": ["6/7=3", "2/9=2000"],
        }
        ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results, as specified in `quacc.schemas.cclib.cclibSchema`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": multiprocessing.cpu_count(),
        "xc": xc,
        "basis": basis,
        "charge": get_charge(atoms) if charge is None else charge,
        "mult": get_multiplicity(atoms) if multiplicity is None else multiplicity,
        "sp": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "pop": "CM5",
        "gfinput": "",
        "ioplist": ["6/7=3", "2/9=2000"],  # see ASE issue #660
    }
    flags = merge_dicts(defaults, calc_swaps)

    atoms.calc = Gaussian(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    return summarize_run(
        atoms,
        LOG_FILE,
        additional_fields={"name": "Gaussian Static"},
    )


@job
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
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
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
        Dictionary of custom kwargs for the calculator. Overrides the following
        defaults:

        ```python
        {
            "mem": "16GB",
            "chk": "Gaussian.chk",
            "nprocshared": multiprocessing.cpu_count(),
            "xc": xc,
            "basis": basis,
            "charge": get_charge(atoms) if charge is None else charge,
            "mult": get_multiplicity(atoms) if multiplicity is None else multiplicity,
            "opt": "",
            "pop": "CM5",
            "scf": ["maxcycle=250", "xqc"],
            "integral": "ultrafine",
            "nosymmetry": "",
            "freq": "" if freq else None,
            "ioplist": ["2/9=2000"],
        }
        ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results, as specified in `quacc.schemas.cclib.cclibSchema`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": multiprocessing.cpu_count(),
        "xc": xc,
        "basis": basis,
        "charge": get_charge(atoms) if charge is None else charge,
        "mult": get_multiplicity(atoms) if multiplicity is None else multiplicity,
        "opt": "",
        "pop": "CM5",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "freq": "" if freq else None,
        "ioplist": ["2/9=2000"],  # ASE issue #660
    }
    flags = merge_dicts(defaults, calc_swaps)

    atoms.calc = Gaussian(**flags)
    atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    return summarize_run(
        atoms,
        LOG_FILE,
        additional_fields={"name": "Gaussian Relax"},
    )
