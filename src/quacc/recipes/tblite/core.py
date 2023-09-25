"""Core recipes for the tblite code"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE
from monty.dev import requires

from quacc import job
from quacc.runners.calc import run_ase_opt, run_ase_vib, run_calc
from quacc.runners.thermo import ideal_gas
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import (
    summarize_opt_run,
    summarize_run,
    summarize_thermo,
    summarize_vib_run,
)
from quacc.utils.dicts import merge_dicts

try:
    from tblite.ase import TBLite
except ImportError:
    TBLite = None

if TYPE_CHECKING:
    from typing import Literal

    from ase import Atoms

    from quacc.schemas.ase import FreqSchema, OptSchema, RunSchema


@job
@requires(TBLite, "tblite must be installed. Refer to the quacc documentation.")
def static_job(
    atoms: Atoms | dict,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    ??? Note

        Calculator Defaults:

        ```python
        {"method": method}
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    method
        GFN1-xTB, GFN2-xTB, and IPEA1-xTB.
    calc_swaps
        Dictionary of custom kwargs for the tblite calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    atoms = fetch_atoms(atoms)

    defaults = {"method": method}
    flags = merge_dicts(defaults, calc_swaps)
    atoms.calc = TBLite(**flags)

    final_atoms = run_calc(atoms, copy_files=copy_files)
    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "TBLite Static"},
    )


@job
@requires(TBLite, "tblite must be installed. Refer to the quacc documentation.")
def relax_job(
    atoms: Atoms | dict,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    relax_cell: bool = False,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Relax a structure.

    ??? Note

        Calculator Defaults:

        ```python
        {"method": method}
        ```

        Optimizer Defaults:

        ```python
        {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB.
    relax_cell
        Whether to relax the cell.
    calc_swaps
        Dictionary of custom kwargs for the tblite calculator.
    opt_swaps
        Dictionary of custom kwargs for [quacc.runners.calc.run_ase_opt][].
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results from [quacc.schemas.ase.summarize_opt_run][]
    """
    atoms = fetch_atoms(atoms)

    defaults = {"method": method}
    flags = merge_dicts(defaults, calc_swaps)
    atoms.calc = TBLite(**flags)

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    dyn = run_ase_opt(atoms, relax_cell=relax_cell, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "TBLite Relax"})


@job
@requires(TBLite, "tblite must be installed. Refer to the quacc documentation.")
def freq_job(
    atoms: Atoms | dict,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    calc_swaps: dict | None = None,
    vib_kwargs: dict | None = None,
    copy_files: list[str] | None = None,
) -> FreqSchema:
    """
    Run a frequency job and calculate thermochemistry.

    ??? Note

        Calculator Defaults:

        ```python
        {"method": method}
        ```

        Frequency Defaults:

        ```python
        {}
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    method
        GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF.
    energy
        Potential energy in eV. If 0, then the output is just the correction.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    calc_swaps
        dictionary of custom kwargs for the tblite calculator.
    vib_kwargs
        dictionary of custom kwargs for the Vibrations object.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    FreqSchema
        Dictionary of results from [quacc.schemas.ase.summarize_vib_run] and
        [quacc.schemas.ase.summarize_thermo][]
    """
    atoms = fetch_atoms(atoms)
    vib_kwargs = vib_kwargs or {}

    defaults = {"method": method}
    flags = merge_dicts(defaults, calc_swaps)
    atoms.calc = TBLite(**flags)

    vibrations = run_ase_vib(atoms, vib_kwargs=vib_kwargs, copy_files=copy_files)
    vib_summary = summarize_vib_run(
        vibrations, additional_fields={"name": "TBLite Frequency"}
    )

    igt = ideal_gas(atoms, vibrations.get_frequencies(), energy=energy)
    vib_summary["thermo"] = summarize_thermo(
        igt,
        temperature=temperature,
        pressure=pressure,
        additional_fields={"name": "ASE Thermo Analysis"},
    )

    return vib_summary
