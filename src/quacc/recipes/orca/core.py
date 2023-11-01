"""Core recipes for ORCA."""
from __future__ import annotations

import multiprocessing
from shutil import which
from typing import TYPE_CHECKING

from ase.calculators.orca import ORCA, OrcaProfile

from quacc import SETTINGS, job
from quacc.runners.calc import run_ase_calc
from quacc.schemas.cclib import cclib_summarize_run
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas.cclib import cclibSchema


LOG_FILE = f"{ORCA().name}.out"
GEOM_FILE = f"{ORCA().name}.xyz"


@job
def static_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    input_swaps: dict[str, Any] | None = None,
    block_swaps: dict[str, Any] | None = None,
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
    input_swaps
        Dictionary of `orcasimpleinput` swaps for the calculator. To enable new
        entries, set the value as True. To remove entries from the defaults, set
        the value as None. For a list of available keys, refer to the
        `ase.calculators.orca.ORCA` calculator.

        !!! Info "Calculator `orcasimpleinput` defaults`"

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
        Dictionary of `orcablock` swaps for the calculator. To enable new entries,
        set the value as True. To remove entries from the defaults, set the
        value as None. For a list of available keys, refer to the
        `ase.calculators.orca.ORCA` calculator.

        !!! Info "Calculator `orcablocks` defaults"

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
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][]
    """

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

    return _base_job(
        atoms,
        charge,
        spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=input_swaps,
        block_swaps=block_swaps,
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
    input_swaps: dict[str, Any] | None = None,
    block_swaps: dict[str, Any] | None = None,
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
    input_swaps
        Dictionary of `orcasimpleinput` swaps for the calculator. To enable new
        entries, set the value as True. To remove entries from the defaults, set
        the value as None. For a list of available keys, refer to the
        `ase.calculators.orca.ORCA` calculator.

        !!! Info "Calculator `orcasimpleinput` defaults"

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
        Dictionary of `orcablock` swaps for the calculator. To enable new entries,
        set the value as True. To remove entries from the defaults, set the
        value as None. For a list of available keys, refer to the
        `ase.calculators.orca.ORCA` calculator.

        !!! Info "Calculator `orcablocks` defaults"

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
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][]
    """

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

    return _base_job(
        atoms,
        charge,
        spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=input_swaps,
        block_swaps=block_swaps,
        additional_fields={"name": "ORCA Relax"},
        copy_files=copy_files,
    )


def _base_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    default_inputs: dict[str, Any] | None = None,
    default_blocks: dict[str, Any] | None = None,
    input_swaps: dict[str, Any] | None = None,
    block_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> cclibSchema:
    """
    Base job function for ORCA recipes.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Charge of the system.
    spin_multiplicity
        Multiplicity of the system.
    default_inputs
        Default input parameters.
    default_blocks
        Default block input parameters.
    input_swaps
        Dictionary of orcasimpleinput swaps for the calculator. To enable new
        entries, set the value as True. To remove entries from the defaults, set
        the value as None.
    block_swaps
        Dictionary of orcablock swaps for the calculator. To enable new entries,
        set the value as True. To remove entries from the defaults, set the
        value as None.
    additional_fields
        Any additional fields to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][]
    """
    inputs = merge_dicts(default_inputs, input_swaps)
    blocks = merge_dicts(default_blocks, block_swaps)
    orcasimpleinput = " ".join(list(inputs.keys()))
    orcablocks = " ".join(list(blocks.keys()))

    atoms.calc = ORCA(
        profile=OrcaProfile([SETTINGS.ORCA_CMD]),
        charge=charge,
        mult=spin_multiplicity,
        orcasimpleinput=orcasimpleinput,
        orcablocks=orcablocks,
    )
    atoms = run_ase_calc(atoms, geom_file=GEOM_FILE, copy_files=copy_files)

    return cclib_summarize_run(
        atoms,
        LOG_FILE,
        additional_fields=additional_fields,
    )
