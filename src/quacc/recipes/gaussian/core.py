"""Core recipes for Gaussian."""
from __future__ import annotations

from typing import TYPE_CHECKING

import psutil

from quacc import job
from quacc.recipes.gaussian._base import base_fn

if TYPE_CHECKING:
    from pathlib import Path

    from ase.atoms import Atoms

    from quacc.schemas._aliases.cclib import cclibSchema


@job
def static_job(
    atoms: Atoms,
    charge: int = 0,
    spin_multiplicity: int = 1,
    xc: str = "wb97xd",
    basis: str = "def2tzvp",
    copy_files: str | Path | list[str | Path] | None = None,
    **calc_kwargs,
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
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the Gaussian calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.gaussian.Gaussian` calculator.

    Returns
    -------
    cclibSchema
        Dictionary of results, as specified in
        [quacc.schemas.cclib.cclib_summarize_run][]
    """

    calc_defaults = {
        "mem": "16GB",
        "chk": "Gaussian.chk",
        "nprocshared": psutil.cpu_count(logical=False),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": spin_multiplicity,
        "sp": "",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
        "pop": "CM5",
        "gfinput": "",
        "ioplist": ["6/7=3", "2/9=2000"],  # see ASE issue #660
    }
    return base_fn(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "Gaussian Static"},
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
    copy_files: str | Path | list[str | Path] | None = None,
    **calc_kwargs,
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
    freq
        If a frequency calculation should be carried out.
    copy_files
        File(s) to copy to the runtime directory. If a directory is provided, it will be recursively unpacked.
    **calc_kwargs
        Custom kwargs for the Gaussian calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.gaussian.Gaussian` calculator.

    Returns
    -------
    cclibSchema
        Dictionary of results, as specified in
        [quacc.schemas.cclib.cclib_summarize_run][]
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
        "freq": "" if freq else None,
        "ioplist": ["2/9=2000"],  # ASE issue #660
    }
    return base_fn(
        atoms,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "Gaussian Relax"},
        copy_files=copy_files,
    )
