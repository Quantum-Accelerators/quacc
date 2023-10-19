"""
Core recipes for Lennard-Jones Potential

NOTE: This set of minimal recipes is mainly for demonstration purposes
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.lj import LennardJones
from ase.optimize import FIRE

from quacc import job
from quacc.builders.thermo import build_ideal_gas
from quacc.runners.calc import run_ase_opt, run_ase_vib, run_calc
from quacc.schemas.ase import summarize_opt_run, summarize_run, summarize_vib_and_thermo
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import OptSchema, RunSchema, VibThermoSchema


@job
def static_job(
    atoms: Atoms,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Function to carry out a static calculation.

    ??? Note

        Calculator Defaults:

        ```python
        {}
        ```

    Parameters
    ----------
    atoms
        Atoms object
    calc_swaps
        Dictionary of custom kwargs for the LJ calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """
    calc_swaps = calc_swaps or {}

    atoms.calc = LennardJones(**calc_swaps)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms, input_atoms=atoms, additional_fields={"name": "LJ Static"}
    )


@job
def relax_job(
    atoms: Atoms,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Function to carry out a geometry optimization.

    ??? Note

        Calculator Defaults:

        ```python
        {}
        ```

        Optimizer Defaults:

        ```python
        {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
        ```

    Parameters
    ----------
    atoms
        Atoms object
    calc_swaps
        Dictionary of custom kwargs for the LJ calculator.
    opt_swaps
        Dictionary of swaps for [quacc.runners.calc.run_ase_opt][].
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """
    calc_swaps = calc_swaps or {}

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    atoms.calc = LennardJones(**calc_swaps)
    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "LJ Relax"})


@job
def freq_job(
    atoms: Atoms,
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    calc_swaps: dict | None = None,
    vib_kwargs: dict | None = None,
    copy_files: list[str] | None = None,
) -> VibThermoSchema:
    """
    Run a frequency job and calculate thermochemistry.

    ??? Note

        Calculator Defaults:

        ```python
        {}
        ```

        Vibrations Defaults:

        ```python
        {}
        ```

    Parameters
    ----------
    atoms
        Atoms object
    energy
        Potential energy in eV. If 0, then the output is just the correction.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    calc_swaps
        dictionary of custom kwargs for the LJ calculator.
    vib_kwargs
        dictionary of custom kwargs for the Vibrations object.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    VibThermoSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_vib_and_thermo][]
    """
    calc_swaps = calc_swaps or {}
    vib_kwargs = vib_kwargs or {}

    atoms.calc = LennardJones(**calc_swaps)
    vibrations = run_ase_vib(atoms, vib_kwargs=vib_kwargs, copy_files=copy_files)
    igt = build_ideal_gas(atoms, vibrations.get_frequencies(), energy=energy)

    return summarize_vib_and_thermo(
        vibrations,
        igt,
        temperature=temperature,
        pressure=pressure,
        additional_fields={"name": "LJ Frequency and Thermo"},
    )
