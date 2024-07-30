"""Core recipes for ORCA."""

from __future__ import annotations

from typing import TYPE_CHECKING

import psutil

from quacc import job
from quacc.atoms.core import perturb
from quacc.recipes.orca._base import run_and_summarize, run_and_summarize_opt

if TYPE_CHECKING:
    from typing import Literal

    from ase.atoms import Atoms
    from numpy.typing import NDArray

    from quacc.types import (
        Filenames,
        OptParams,
        SourceDirectory,
        cclibASEOptSchema,
        cclibSchema,
    )


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    nprocs: int | Literal["max"] = "max",
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
        Dictionary of results from [quacc.schemas.cclib.CclibSummarize.run][].
        See the type-hint for the data structure.
    """
    nprocs = psutil.cpu_count(logical=False) if nprocs == "max" else nprocs
    default_inputs = [xc, basis, "engrad", "normalprint"]
    default_blocks = [f"%pal nprocs {nprocs} end"]

    return run_and_summarize(
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
    nprocs: int | Literal["max"] = "max",
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
        Dictionary of results from [quacc.schemas.cclib.CclibSummarize.run][].
        See the type-hint for the data structure.
    """
    nprocs = psutil.cpu_count(logical=False) if nprocs == "max" else nprocs

    default_inputs = [xc, basis, "normalprint", "opt"]
    if run_freq:
        default_inputs.append("freq")

    default_blocks = [f"%pal nprocs {nprocs} end"]

    return run_and_summarize(
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
def freq_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    numerical: bool = False,
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    nprocs: int | Literal["max"] = "max",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> cclibSchema:
    """
    Carry out a vibrational frequency analysis calculation.

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
    numerical
        If True (default False), a numeric frequency calculation will be requested
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
        Dictionary of results from [quacc.schemas.cclib.CclibSummarize.run][].
        See the type-hint for the data structure.
    """
    nprocs = psutil.cpu_count(logical=False) if nprocs == "max" else nprocs

    default_inputs = [xc, basis, "normalprint", "numfreq" if numerical else "freq"]

    default_blocks = [f"%pal nprocs {nprocs} end"]

    return run_and_summarize(
        atoms,
        charge,
        spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=orcasimpleinput,
        block_swaps=orcablocks,
        additional_fields={"name": "ORCA vibrational frequency analysis"},
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
    opt_params: OptParams | None = None,
    nprocs: int | Literal["max"] = "max",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> cclibASEOptSchema:
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
        Dictionary of results from [quacc.schemas.cclib.CclibSummarize.run][] merged with
        the results from [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    nprocs = psutil.cpu_count(logical=False) if nprocs == "max" else nprocs
    default_inputs = [xc, basis, "engrad", "normalprint"]
    default_blocks = [f"%pal nprocs {nprocs} end"]

    return run_and_summarize_opt(
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


@job
def ase_quasi_irc_job(
    atoms: Atoms,
    mode: list[list[float]] | NDArray,
    perturb_magnitude: float = 0.6,
    direction: Literal["forward", "reverse"] = "forward",
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97x-d3bj",
    basis: str = "def2-tzvp",
    orcasimpleinput: list[str] | None = None,
    orcablocks: list[str] | None = None,
    opt_params: OptParams | None = None,
    nprocs: int | Literal["max"] = "max",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> cclibASEOptSchema:
    """
    Quasi-IRC to optimize a reaction endpoint from a transition-state with known vibrational frequency modes.
    Perturbs the structure of `atoms` by a finite amount (0.6 * the normalized mode magnitude) along the specified
    vibrational frequency mode (assumed to be the transition mode), and then performs a `relax_job` on the perturbed
    structure.

    Parameters
    ----------
    atoms
        Atoms object
    mode
        Transition mode. This should be an Nx3 matrix, where N is the number of atoms in `atoms`.
    perturb_magnitude
        Factor to multiply the transition mode. Default is 0.6. In some cases, it may be advisable to increase this
        factor, perhaps to 1.0 or 1.1. Lowering it is not generally found to be helpful.
    direction
        Direction of the (Quasi)IRC. Should be "forward" or "reverse".
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
        Dictionary of results from [quacc.schemas.cclib.CclibSummarize.run][] merged with
        the results from [quacc.schemas.ase.Summarize.opt][].
        See the type-hint for the data structure.
    """
    nprocs = psutil.cpu_count(logical=False) if nprocs == "max" else nprocs
    default_inputs = [xc, basis, "engrad", "normalprint"]
    default_blocks = [f"%pal nprocs {nprocs} end"]

    scale = perturb_magnitude if direction == "forward" else perturb_magnitude * -1

    return run_and_summarize_opt(
        perturb(atoms, mode, scale),
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        default_inputs=default_inputs,
        default_blocks=default_blocks,
        input_swaps=orcasimpleinput,
        block_swaps=orcablocks,
        opt_params=opt_params,
        additional_fields={"name": "ORCA ASE Quasi-IRC optimization"},
        copy_files=copy_files,
    )
