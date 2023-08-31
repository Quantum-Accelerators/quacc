"""
Core recipes for Lennard-Jones Potential

NOTE: This set of minimal recipes is mainly for demonstration purposes
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.lj import LennardJones
from ase.optimize import FIRE

from quacc import job
from quacc.schemas.ase import (
    summarize_opt_run,
    summarize_run,
    summarize_thermo,
    summarize_vib_run,
)
from quacc.utils.calc import run_ase_opt, run_ase_vib, run_calc
from quacc.utils.dicts import merge_dicts
from quacc.utils.thermo import ideal_gas
from quacc.utils.wflows import fetch_atoms

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import OptSchema, RunSchema, ThermoSchema, VibSchema

    class FreqSchema(VibSchema):
        thermo: ThermoSchema


@job
def static_job(
    atoms: Atoms | dict,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Function to carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    calc_swaps
        Dictionary of custom kwargs for the LJ calculator
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in `quacc.schemas.ase.RunSchema`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    atoms.calc = LennardJones(**calc_swaps)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms, input_atoms=atoms, additional_fields={"name": "LJ Static"}
    )


@job
def relax_job(
    atoms: Atoms | dict,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Function to carry out a geometry optimization

    Parameters
    ----------
    atoms
        Atoms object
    calc_swaps
        Dictionary of custom kwargs for the LJ calculator.
    opt_swaps
        Dictionary of swaps for run_ase_opt
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in `quacc.schemas.ase.OptSchema`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}

    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    atoms.calc = LennardJones(**calc_swaps)
    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "LJ Relax"})


@job
def freq_job(
    atoms: Atoms | dict,
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
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    energy
        Potential energy in eV. If 0, then the output is just the correction.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    calc_swaps
        dictionary of custom kwargs for the LJ calculator.
    vib_kwargs
        dictionary of custom kwargs for the Vibrations object
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    FreqSchema
        Dictionary of results, specified in `quacc.schemas.ase.FreqSchema`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}
    vib_kwargs = vib_kwargs or {}

    atoms.calc = LennardJones(**calc_swaps)
    vibrations = run_ase_vib(atoms, vib_kwargs=vib_kwargs, copy_files=copy_files)
    vib_summary = summarize_vib_run(
        vibrations, additional_fields={"name": "LJ Frequency Analysis"}
    )

    igt = ideal_gas(atoms, vibrations.get_frequencies(), energy=energy)
    vib_summary["thermo"] = summarize_thermo(
        igt, temperature=temperature, pressure=pressure
    )

    return vib_summary
