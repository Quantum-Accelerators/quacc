"""Core recipes for ORCA"""
from __future__ import annotations

import multiprocessing
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.orca import ORCA
import covalent as ct
from quacc.schemas.cclib import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc

LOG_FILE = ORCA().name + ".out"
GEOM_FILE = ORCA().name + ".xyz"


@ct.electron
def StaticJob(
    atoms: Atoms,
    charge: int = None,
    mult: int = None,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    input_swaps: Dict[str, Any] | None = None,
    block_swaps: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
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
    input_swaps
        Dictionary of orcasimpleinput swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.
    block_swaps
        Dictionary of orcablock swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.

    Returns
    -------
    summary
        Summary of the run.
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

    inputs = merge_dicts(
        default_inputs, input_swaps, remove_none=True, remove_false=True
    )
    blocks = merge_dicts(
        default_blocks, block_swaps, remove_none=True, remove_false=True
    )
    orcasimpleinput = " ".join(list(inputs.keys()))
    orcablocks = " ".join(list(blocks.keys()))

    atoms.calc = ORCA(
        charge=charge if charge else round(sum(atoms.get_initial_charges())),
        mult=mult if mult else round(1 + sum(atoms.get_initial_magnetic_moments())),
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_calc(atoms, geom_file=GEOM_FILE)
    summary = summarize_run(atoms, LOG_FILE)

    return summary


@ct.electron
def RelaxJob(
    atoms: Atoms,
    charge: int = None,
    mult: int = None,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    freq: bool = False,
    input_swaps: Dict[str, Any] | None = None,
    block_swaps: Dict[str, Any] | None = None,
):
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
    input_swaps
        Dictionary of orcasimpleinput swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.
    block_swaps
        Dictionary of orcablock swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None/False.

    Returns
    -------
    summary
        Summary of the run.
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
        "freq": True if freq else None,
        "xyzfile": True,
    }
    default_blocks = {}

    inputs = merge_dicts(
        default_inputs, input_swaps, remove_none=True, remove_false=True
    )
    blocks = merge_dicts(
        default_blocks, block_swaps, remove_none=True, remove_false=True
    )
    orcasimpleinput = " ".join(list(inputs.keys()))
    orcablocks = " ".join(list(blocks.keys()))

    atoms.calc = ORCA(
        charge=charge if charge else round(sum(atoms.get_initial_charges())),
        mult=mult if mult else round(1 + sum(atoms.get_initial_magnetic_moments())),
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_calc(atoms, geom_file=GEOM_FILE)
    summary = summarize_run(atoms, LOG_FILE)

    return summary
