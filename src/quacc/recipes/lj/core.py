"""
Core recipes for Lennard-Jones Potential.

NOTE: This set of minimal recipes is mainly for demonstration purposes
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.lj import LennardJones
from ase.optimize import FIRE

from quacc import job
from quacc.runners.ase import run_calc, run_opt, run_vib
from quacc.runners.thermo import run_ideal_gas
from quacc.schemas.ase import summarize_opt_run, summarize_run, summarize_vib_and_thermo
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.runners.ase import VibKwargs
    from quacc.schemas._aliases.ase import OptSchema, RunSchema, VibThermoSchema


@job
def static_job(atoms: Atoms, **calc_kwargs) -> RunSchema:
    """
    Function to carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    **calc_kwargs
        Dictionary of custom kwargs for the LJ calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.lj.LJ` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    atoms.calc = LennardJones(**calc_kwargs)
    final_atoms = run_calc(atoms)

    return summarize_run(final_atoms, atoms, additional_fields={"name": "LJ Static"})


@job
def relax_job(
    atoms: Atoms, opt_params: dict[str, Any] | None = None, **calc_kwargs
) -> OptSchema:
    """
    Function to carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    **calc_kwargs
        Custom kwargs for the LJ calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.lj.LJ` calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """
    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": FIRE}
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    atoms.calc = LennardJones(**calc_kwargs)
    dyn = run_opt(atoms, **opt_flags)

    return summarize_opt_run(dyn, additional_fields={"name": "LJ Relax"})


@job
def freq_job(
    atoms: Atoms,
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    vib_kwargs: VibKwargs | None = None,
    **calc_kwargs,
) -> VibThermoSchema:
    """
    Run a frequency job and calculate thermochemistry.

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
    vib_kwargs
        Dictionary of kwargs for the `ase.vibrations.Vibrations` class.
    **calc_kwargs
        Dictionary of custom kwargs for the LJ calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.lj.LJ` calculator.

    Returns
    -------
    VibThermoSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_vib_and_thermo][].
        See the type-hint for the data structure.
    """
    vib_kwargs = vib_kwargs or {}

    atoms.calc = LennardJones(**calc_kwargs)
    vibrations = run_vib(atoms, vib_kwargs=vib_kwargs)
    igt = run_ideal_gas(atoms, vibrations.get_frequencies(), energy=energy)

    return summarize_vib_and_thermo(
        vibrations,
        igt,
        temperature=temperature,
        pressure=pressure,
        additional_fields={"name": "LJ Frequency and Thermo"},
    )
