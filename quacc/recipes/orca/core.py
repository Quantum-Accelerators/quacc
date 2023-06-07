"""Core recipes for ORCA"""
from __future__ import annotations

import multiprocessing

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.orca import ORCA

from quacc.schemas.cclib import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import remove_dict_empties

LOG_FILE = f"{ORCA().name}.out"
GEOM_FILE = f"{ORCA().name}.xyz"


@ct.electron
def static_job(
    atoms: Atoms,
    charge: int | None = None,
    mult: int | None = None,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    input_swaps: dict | None = None,
    block_swaps: dict | None = None,
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
    input_swaps
        dictionary of orcasimpleinput swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.
    block_swaps
        dictionary of orcablock swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.cclib.summarize_run
    """

    input_swaps = input_swaps or {}
    block_swaps = block_swaps or {}

    if not any(k for k in block_swaps if "nprocs" in k.lower()):
        nprocs = multiprocessing.cpu_count()
        block_swaps[f"%pal nprocs {nprocs} end"] = True

    default_inputs = {
        xc: True,
        basis: True,
        "sp": True,
        "slowconv": True,
        "normalprint": True,
        "xyzfile": True,
    }
    default_blocks = {}

    inputs = remove_dict_empties(default_inputs | input_swaps)
    blocks = remove_dict_empties(default_blocks | block_swaps)
    orcasimpleinput = " ".join(list(inputs.keys()))
    orcablocks = " ".join(list(blocks.keys()))

    atoms.calc = ORCA(
        charge=charge or round(sum(atoms.get_initial_charges())),
        mult=mult or round(1 + sum(atoms.get_initial_magnetic_moments())),
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_calc(atoms, geom_file=GEOM_FILE)

    return summarize_run(atoms, LOG_FILE, additional_fields={"name": "ORCA Static"})


@ct.electron
def relax_job(
    atoms: Atoms,
    charge: int = None,
    mult: int = None,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    run_freq: bool = False,
    input_swaps: dict = None,
    block_swaps: dict = None,
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
    run_freq
        If a requency calculation should be carried out.
    input_swaps
        dictionary of orcasimpleinput swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.
    block_swaps
        dictionary of orcablock swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.cclib.summarize_run
    """

    input_swaps = input_swaps or {}
    block_swaps = block_swaps or {}

    if not any(k for k in block_swaps if "nprocs" in k.lower()):
        nprocs = multiprocessing.cpu_count()
        block_swaps[f"%pal nprocs {nprocs} end"] = True

    default_inputs = {
        xc: True,
        basis: True,
        "opt": True,
        "slowconv": True,
        "normalprint": True,
        "freq": True if run_freq else None,
        "xyzfile": True,
    }
    default_blocks = {}

    inputs = remove_dict_empties(default_inputs | input_swaps)
    blocks = remove_dict_empties(default_blocks | block_swaps)
    orcasimpleinput = " ".join(list(inputs.keys()))
    orcablocks = " ".join(list(blocks.keys()))

    atoms.calc = ORCA(
        charge=charge or round(sum(atoms.get_initial_charges())),
        mult=mult or round(1 + sum(atoms.get_initial_magnetic_moments())),
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_calc(atoms, geom_file=GEOM_FILE)

    return summarize_run(atoms, LOG_FILE, additional_fields={"name": "ORCA Relax"})
