"""Core recipes for ORCA"""
from __future__ import annotations

import multiprocessing
import os
from typing import TYPE_CHECKING

import covalent as ct
from ase.calculators.orca import ORCA, OrcaProfile

from quacc import SETTINGS
from quacc.schemas.atoms import fetch_atoms
from quacc.schemas.cclib import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import get_parameters

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.cclib import cclibSchema


LOG_FILE = f"{ORCA().name}.out"
GEOM_FILE = f"{ORCA().name}.xyz"


@ct.electron
def static_job(
    atoms: Atoms | dict,
    charge: int | None = None,
    multiplicity: int | None = None,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    input_swaps: dict | None = None,
    block_swaps: dict | None = None,
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
        `atoms.get_initial_charges()`.
    multiplicity
        Multiplicity of the system. If None, this is determined from 1+ the sum
        of `atoms.get_initial_magnetic_moments()`.
    xc
        Exchange-correlation functional
    basis
        Basis set
    input_swaps
        Dictionary of orcasimpleinput swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None.
    block_swaps
        Dictionary of orcablock swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results from quacc.schemas.cclib.summarize_run
    """

    atoms = fetch_atoms(atoms)
    block_swaps = block_swaps or {}

    default_inputs = {
        xc: True,
        basis: True,
        "sp": True,
        "slowconv": True,
        "normalprint": True,
        "xyzfile": True,
    }
    default_blocks = {}

    if not any(k for k in block_swaps if "nprocs" in k.lower()) and os.environ.get(
        "mpirun"
    ):
        nprocs = multiprocessing.cpu_count()
        block_swaps[f"%pal nprocs {nprocs} end"] = True

    inputs = get_parameters(default_inputs, swaps=input_swaps)
    blocks = get_parameters(default_blocks, swaps=block_swaps)

    orcasimpleinput = " ".join(list(inputs.keys()))
    orcablocks = " ".join(list(blocks.keys()))

    charge = int(atoms.get_initial_charges().sum()) if charge is None else charge
    multiplicity = (
        int(1 + atoms.get_initial_magnetic_moments().sum())
        if multiplicity is None
        else multiplicity
    )

    calc = ORCA(
        profile=OrcaProfile([SETTINGS.ORCA_CMD]),
        charge=charge,
        mult=multiplicity,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_calc(atoms, calc, geom_file=GEOM_FILE, copy_files=copy_files)

    return summarize_run(
        atoms,
        LOG_FILE,
        additional_fields={"name": "ORCA Static"},
    )


@ct.electron
def relax_job(
    atoms: Atoms | dict,
    charge: int | None = None,
    multiplicity: int | None = None,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    run_freq: bool = False,
    input_swaps: dict | None = None,
    block_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> cclibSchema:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system. If None, this is determined from the sum of
        atoms.get_initial_charges().
    multiplicity
        Multiplicity of the system. If None, this is determined from 1+ the sum
        of atoms.get_initial_magnetic_moments().
    xc
        Exchange-correlation functional
    basis
        Basis set
    run_freq
        If a requency calculation should be carried out.
    input_swaps
        Dictionary of orcasimpleinput swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None.
    block_swaps
        Dictionary of orcablock swaps for the calculator.
        To enable new entries, set the value as True.
        To remove entries from the defaults, set the value as None.
    copy_files
        Absolute paths to files to copy to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results from quacc.schemas.cclib.summarize_run
    """

    atoms = fetch_atoms(atoms)
    block_swaps = block_swaps or {}

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

    if not any(k for k in block_swaps if "nprocs" in k.lower()) and os.environ.get(
        "mpirun"
    ):
        nprocs = multiprocessing.cpu_count()
        block_swaps[f"%pal nprocs {nprocs} end"] = True

    inputs = get_parameters(default_inputs, swaps=input_swaps)
    blocks = get_parameters(default_blocks, swaps=block_swaps)

    orcasimpleinput = " ".join(list(inputs.keys()))
    orcablocks = " ".join(list(blocks.keys()))

    charge = int(atoms.get_initial_charges().sum()) if charge is None else charge
    multiplicity = (
        int(1 + atoms.get_initial_magnetic_moments().sum())
        if multiplicity is None
        else multiplicity
    )

    calc = ORCA(
        profile=OrcaProfile([SETTINGS.ORCA_CMD]),
        charge=charge,
        mult=multiplicity,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_calc(atoms, calc, geom_file=GEOM_FILE, copy_files=copy_files)

    return summarize_run(
        atoms,
        LOG_FILE,
        additional_fields={"name": "ORCA Relax"},
    )
