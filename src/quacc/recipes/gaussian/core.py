"""Core recipes for Gaussian."""

from __future__ import annotations

from typing import TYPE_CHECKING

import psutil

from quacc import job
from quacc.recipes.gaussian._base import run_and_summarize

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, RunSchema, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97xd",
    basis: str = "def2tzvp",
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
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
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the Gaussian calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.gaussian.Gaussian][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results
    """
    calc_defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": psutil.cpu_count(logical=False),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": spin_multiplicity,
        "force": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "pop": "CM5",
        "gfinput": "",
        "ioplist": ["6/7=3", "2/9=2000"],  # see ASE issue #660
    }
    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "Gaussian Static"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    charge: int,
    spin_multiplicity: int,
    xc: str = "wb97xd",
    basis: str = "def2tzvp",
    freq: bool = False,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
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
    freq
        If a frequency calculation should be carried out.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    additional_fields
        Additional fields to add to the results
    **calc_kwargs
        Custom kwargs for the Gaussian calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.gaussian.Gaussian][] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results
    """
    calc_defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": psutil.cpu_count(logical=False),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": spin_multiplicity,
        "opt": "",
        "pop": "CM5",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "ioplist": ["2/9=2000"],  # ASE issue #660
    }
    if freq:
        calc_defaults["freq"] = ""

    return run_and_summarize(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "Gaussian Relax"} | (additional_fields or {}),
        copy_files=copy_files,
    )
