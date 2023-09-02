"""Core recipes for ORCA"""
from __future__ import annotations

import multiprocessing
from shutil import which
from typing import TYPE_CHECKING

from ase.calculators.orca import ORCA, OrcaProfile

from quacc import SETTINGS, job
from quacc.schemas import fetch_atoms
from quacc.schemas.cclib import summarize_run
from quacc.utils.atoms import get_charge, get_multiplicity
from quacc.utils.calc import run_calc
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.cclib import cclibSchema


LOG_FILE = f"{ORCA().name}.out"
GEOM_FILE = f"{ORCA().name}.xyz"


@job
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
    input_swaps
        Dictionary of orcasimpleinput swaps for the calculator. To enable new
        entries, set the value as True. To remove entries from the defaults, set
        the value as None. Overrides the following defaults:

        ```python
        {
            xc: True,
            basis: True,
            "sp": True,
            "slowconv": True,
            "normalprint": True,
            "xyzfile": True,
        }
        ```
    block_swaps
        Dictionary of orcablock swaps for the calculator. To enable new entries,
        set the value as True. To remove entries from the defaults, set the
        value as None. Overrides the following defaults:

        ```python
        (
            {f"%pal nprocs {multiprocessing.cpu_count()} end": True}
            if which("mpirun")
            else {}
        )
        ```

    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results from `quacc.schemas.cclib.summarize_run`
    """
    atoms = fetch_atoms(atoms)
    input_swaps = input_swaps or {}
    block_swaps = block_swaps or {}

    default_inputs = {
        xc: True,
        basis: True,
        "sp": True,
        "slowconv": True,
        "normalprint": True,
        "xyzfile": True,
    }
    default_blocks = (
        {f"%pal nprocs {multiprocessing.cpu_count()} end": True}
        if which("mpirun")
        else {}
    )

    inputs = merge_dicts(default_inputs, input_swaps)
    blocks = merge_dicts(default_blocks, block_swaps)
    orcasimpleinput = " ".join(list(inputs.keys()))
    orcablocks = " ".join(list(blocks.keys()))

    atoms.calc = ORCA(
        profile=OrcaProfile([SETTINGS.ORCA_CMD]),
        charge=get_charge(atoms) if charge is None else charge,
        mult=get_multiplicity(atoms) if multiplicity is None else multiplicity,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    return summarize_run(
        atoms,
        LOG_FILE,
        additional_fields={"name": "ORCA Static"},
    )


@job
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
        Dictionary of orcasimpleinput swaps for the calculator. To enable new
        entries, set the value as True. To remove entries from the defaults, set
        the value as None. Overrides the following defaults:

        ```python
        {
            xc: True,
            basis: True,
            "opt": True,
            "slowconv": True,
            "normalprint": True,
            "freq": True if run_freq else None,
            "xyzfile": True,
        }
        ```
    block_swaps
        Dictionary of orcablock swaps for the calculator. To enable new entries,
        set the value as True. To remove entries from the defaults, set the
        value as None. Overrides the following defaults:

        ```python
        (
            {f"%pal nprocs {multiprocessing.cpu_count()} end": True}
            if which("mpirun")
            else {}
        )
        ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results from `quacc.schemas.cclib.summarize_run`
    """
    atoms = fetch_atoms(atoms)
    input_swaps = input_swaps or {}
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
    default_blocks = (
        {f"%pal nprocs {multiprocessing.cpu_count()} end": True}
        if which("mpirun")
        else {}
    )

    inputs = merge_dicts(default_inputs, input_swaps)
    blocks = merge_dicts(default_blocks, block_swaps)
    orcasimpleinput = " ".join(list(inputs.keys()))
    orcablocks = " ".join(list(blocks.keys()))

    atoms.calc = ORCA(
        profile=OrcaProfile([SETTINGS.ORCA_CMD]),
        charge=get_charge(atoms) if charge is None else charge,
        mult=get_multiplicity(atoms) if multiplicity is None else multiplicity,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    return summarize_run(
        atoms,
        LOG_FILE,
        additional_fields={"name": "ORCA Relax"},
    )
