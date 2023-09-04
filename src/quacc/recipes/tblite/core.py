"""Core recipes for the tblite code"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE
from monty.dev import requires

from quacc import job
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import (
    summarize_opt_run,
    summarize_run,
    summarize_thermo,
    summarize_vib_run,
)
from quacc.utils.calc import run_ase_opt, run_ase_vib, run_calc
from quacc.utils.dicts import merge_dicts
from quacc.utils.thermo import ideal_gas

try:
    from tblite.ase import TBLite
except ImportError:
    TBLite = None

if TYPE_CHECKING:
    from typing import Literal

    from ase import Atoms

    from quacc.schemas.ase import OptSchema, RunSchema, ThermoSchema, VibSchema

    class FreqSchema(VibSchema):
        thermo: ThermoSchema


@job
@requires(TBLite, "tblite must be installed. Try pip install tblite[ase]")
def static_job(
    atoms: Atoms | dict,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    method
        GFN1-xTB, GFN2-xTB, and IPEA1-xTB.
    calc_swaps
        Dictionary of custom kwargs for the tblite calculator. Overrides the
        following defaults: `{}`.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    atoms.calc = TBLite(method=method, **calc_swaps)
    final_atoms = run_calc(atoms, copy_files=copy_files)
    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "TBLite Static"},
    )


@job
@requires(TBLite, "tblite must be installed. Try pip install tblite[ase]")
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
        Dictionary of custom kwargs for the tblite calculator. Overrides the
        following defaults: `{}`
    opt_swaps
        Dictionary of custom kwargs for `run_ase_opt`. Overrides the following
        defaults:

        ```python
        {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
        ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results from `quacc.schemas.ase.summarize_opt_run`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    atoms.calc = TBLite(method=method, **calc_swaps)
    dyn = run_ase_opt(atoms, relax_cell=relax_cell, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "TBLite Relax"})


@job
@requires(TBLite, "tblite must be installed. Try pip install tblite[ase]")
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
        dictionary of custom kwargs for the xTB calculator. Overrides the
        following defaults: `{}`
    vib_kwargs
        dictionary of custom kwargs for the Vibrations object.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    FreqSchema
        Dictionary of results from `quacc.schemas.ase.summarize_vib_run` patched
        with the results of `quacc.schemas.ase.summarize_thermo` in the "thermo"
        key.
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}
    vib_kwargs = vib_kwargs or {}

    atoms.calc = TBLite(method=method, **calc_swaps)
    vibrations = run_ase_vib(atoms, vib_kwargs=vib_kwargs, copy_files=copy_files)
    vib_summary = summarize_vib_run(
        vibrations, additional_fields={"name": "TBLite Frequency Analysis"}
    )

    igt = ideal_gas(atoms, vibrations.get_frequencies(), energy=energy)
    vib_summary["thermo"] = summarize_thermo(
        igt, temperature=temperature, pressure=pressure
    )

    return vib_summary
