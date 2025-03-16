"""
Core recipes for Lennard-Jones Potential.

NOTE: This set of minimal recipes is mainly for demonstration purposes
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.lj import LennardJones
from ase.optimize import BFGS

from quacc import job
from quacc.recipes.common.core import Recipe
from quacc.schemas.thermo import ThermoSummarize
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from ase.optimize.optimize import Optimizer

    from quacc.types import OptSchema, RunSchema, VibKwargs, VibThermoSchema


@job
def static_job(
    atoms: Atoms, additional_fields: dict[str, Any] | None = None, **calc_kwargs
) -> RunSchema:
    """
    Function to carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    additional_fields
        Additional fields for results
    **calc_kwargs
        Dictionary of custom kwargs for the LJ calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.lj.LennardJones] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    return Recipe(LennardJones).run(
        atoms, additional_fields=additional_fields, **calc_kwargs
    )


@job
def relax_job(
    atoms: Atoms,
    fmax: float | None = 0.01,
    max_steps: int = 1000,
    optimizer: type[Optimizer] = BFGS,
    optimizer_kwargs: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Function to carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    fmax
        Maximum force change in eV/A
    max_steps
        Maximum number of steps
    optimizer
        ASE optimizer class to use
    optimizer_kwargs
        Dictionary of keyword arguments to pass to the optimizer
    additional_fields
        Metadata to store in the results
    **calc_kwargs
        Calculator parameters to pass to [ase.calculators.lj.LennardJones][]

    Returns
    -------
    OptSchema
        Results dictionary
    """
    return Recipe(LennardJones).relax(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        optimizer_kwargs=optimizer_kwargs,
        additional_fields=additional_fields,
        **calc_kwargs,
    )


@job
def freq_job(
    atoms: Atoms,
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    vib_kwargs: VibKwargs | None = None,
    additional_fields: dict[str, Any] | None = None,
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
        Dictionary of kwargs for the [ase.vibrations.Vibrations][] class.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Dictionary of custom kwargs for the LJ calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.lj.LennardJones] calculator.

    Returns
    -------
    VibThermoSchema
        Dictionary of results
    """
    freq_results = Recipe(LennardJones).vib(
        atoms,
        energy=energy,
        temperature=temperature,
        pressure=pressure,
        vib_kwargs=vib_kwargs,
        additional_fields=additional_fields,
        **calc_kwargs,
    )
    thermo_results = ThermoSummarize(
        atoms,
        freq_results["results"]["vib_freqs"],
        energy=energy,
        additional_fields=additional_fields,
    ).ideal_gas(temperature=temperature, pressure=pressure)
    return recursive_dict_merge(freq_results, thermo_results)
