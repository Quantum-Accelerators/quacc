"""Core recipes for ORCA."""
from __future__ import annotations

from typing import TYPE_CHECKING

import psutil

from quacc import job
from quacc.recipes.orca._base import base_fn

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.cclib import cclibSchema


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    orcasimpleinput: dict[str, Any] | None = None,
    orcablocks: dict[str, Any] | None = None,
    nprocs: int | None = None,
    copy_files: list[str] | None = None,
) -> cclibSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    xc
        Exchange-correlation functional
    basis
        Basis set
    orcasimpleinput
        Dictionary of `orcasimpleinput` swaps for the calculator. To enable new
        entries, set the value as True. To remove entries from the defaults, set
        the value as None. For a list of available keys, refer to the
        `ase.calculators.orca.ORCA` calculator.
    orcablocks
        Dictionary of `orcablocks` swaps for the calculator. To enable new entries,
        set the value as True. To remove entries from the defaults, set the
        value as None. For a list of available keys, refer to the
        `ase.calculators.orca.ORCA` calculator.
    nprocs
        Number of processors to use. Defaults to the number of physical cores.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][]
    """

    nprocs = nprocs or psutil.cpu_count(logical=False)
    default_inputs = {
        xc: True,
        basis: True,
        "sp": True,
        "slowconv": True,
        "normalprint": True,
        "xyzfile": True,
    }
    default_blocks = {f"%pal nprocs {nprocs} end": True}

    return base_fn(
        atoms,
        charge,
        spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=orcasimpleinput,
        block_swaps=orcablocks,
        additional_fields={"name": "ORCA Static"},
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    run_freq: bool = False,
    orcasimpleinput: dict[str, Any] | None = None,
    orcablocks: dict[str, Any] | None = None,
    nprocs: int | None = None,
    copy_files: list[str] | None = None,
) -> cclibSchema:
    """
    Carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    xc
        Exchange-correlation functional
    basis
        Basis set
    run_freq
        If a frequency calculation should be carried out.
    orcasimpleinput
        Dictionary of `orcasimpleinput` swaps for the calculator. To enable new
        entries, set the value as True. To remove entries from the defaults, set
        the value as None. For a list of available keys, refer to the
        `ase.calculators.orca.ORCA` calculator.
    orcablocks
        Dictionary of `orcablocks` swaps for the calculator. To enable new entries,
        set the value as True. To remove entries from the defaults, set the
        value as None. For a list of available keys, refer to the
        `ase.calculators.orca.ORCA` calculator.
    nprocs
        Number of processors to use. Defaults to the number of physical cores.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][]
    """

    nprocs = nprocs or psutil.cpu_count(logical=False)
    default_inputs = {
        xc: True,
        basis: True,
        "opt": True,
        "slowconv": True,
        "normalprint": True,
        "freq": True if run_freq else None,
        "xyzfile": True,
    }
    default_blocks = {f"%pal nprocs {nprocs} end": True}

    return base_fn(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=orcasimpleinput,
        block_swaps=orcablocks,
        additional_fields={"name": "ORCA Relax"},
        copy_files=copy_files,
    )
