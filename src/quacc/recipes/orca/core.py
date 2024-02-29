"""Core recipes for ORCA."""

from __future__ import annotations

from typing import TYPE_CHECKING

import psutil

from quacc import job
from quacc.recipes.orca._base import base_fn, base_opt_fn

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.cclib import cclibSchema
    from quacc.utils.files import Filenames, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    nprocs: int | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        List of `orcasimpleinput` swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name. Refer to the
        [ase.calculators.orca.ORCA][] calculator for details on `orcasimpleinput`.
    orcablocks
        List of `orcablocks` swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name. Refer to the
        [ase.calculators.orca.ORCA][] calculator for details on `orcablocks`.
    nprocs
        Number of processors to use. Defaults to the number of physical cores.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][].
        See the type-hint for the data structure.
    """

    nprocs = nprocs or psutil.cpu_count(logical=False)
    default_inputs = [xc, basis, "sp", "slowconv", "normalprint", "xyzfile"]
    default_blocks = [f"%pal nprocs {nprocs} end"]

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
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    run_freq: bool = False,
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    nprocs: int | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        List of `orcasimpleinput` swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name. Refer to the
        [ase.calculators.orca.ORCA][] calculator for details on `orcasimpleinput`.
    orcablocks
        List of `orcablocks` swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name. Refer to the
        [ase.calculators.orca.ORCA][] calculator for details on `orcablocks`.
    nprocs
        Number of processors to use. Defaults to the number of physical cores.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    cclibSchema
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][].
        See the type-hint for the data structure.
    """
    nprocs = nprocs or psutil.cpu_count(logical=False)

    default_inputs = [xc, basis, "opt", "slowconv", "normalprint", "xyzfile"]
    if run_freq:
        default_inputs.append("freq")

    default_blocks = [f"%pal nprocs {nprocs} end"]

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


@job
def ase_relax_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    opt_params: dict[str, Any] | None = None,
    nprocs: int | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        Basis set.
    orcasimpleinput
        List of `orcasimpleinput` swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name. Refer to the
        [ase.calculators.orca.ORCA][] calculator for details on `orcasimpleinput`.
    orcablocks
        List of `orcablocks` swaps for the calculator. To remove entries
        from the defaults, put a `#` in front of the name. Refer to the
        [ase.calculators.orca.ORCA][] calculator for details on `orcablocks`.
    nprocs
        Number of processors to use. Defaults to the number of physical cores.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    cclibASEOptSchema
        Dictionary of results from [quacc.schemas.cclib.cclib_summarize_run][] merged with
        the results from [quacc.schemas.ase.summarize_opt_run][].
        See the type-hint for the data structure.
    """

    nprocs = nprocs or psutil.cpu_count(logical=False)
    default_inputs = [xc, basis, "slowconv", "normalprint", "xyzfile", "engrad"]
    default_blocks = [f"%pal nprocs {nprocs} end"]

    return base_opt_fn(
        atoms,
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=orcasimpleinput,
        block_swaps=orcablocks,
        opt_params=opt_params,
        additional_fields={"name": "ORCA ASE Relax"},
        copy_files=copy_files,
    )
