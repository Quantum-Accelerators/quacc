"""
Core recipes for Lennard-Jones Potential.

NOTE: This set of minimal recipes is mainly for demonstration purposes
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.lj import LennardJones

from quacc import job
from quacc.recipes.common.core import Recipe
from quacc.runners.ase import Runner
from quacc.schemas.ase import VibSummarize

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import OptParams, OptSchema, RunSchema, VibKwargs, VibThermoSchema


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
    opt_params: OptParams | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Function to carry out a geometry optimization.

    Parameters
    ----------
    atoms
        Atoms object
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Custom kwargs for the LJ calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.lj.LennardJones] calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    return Recipe(LennardJones).relax(
        atoms, opt_params=opt_params, additional_fields=additional_fields, **calc_kwargs
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
    vib_kwargs = vib_kwargs or {}

    calc = LennardJones(**calc_kwargs)
    vib = Runner(atoms, calc).run_vib(vib_kwargs=vib_kwargs)

    return VibSummarize(
        vib,
        additional_fields={"name": "LJ Frequency and Thermo"}
        | (additional_fields or {}),
    ).vib_and_thermo(
        "ideal_gas", energy=energy, temperature=temperature, pressure=pressure
    )
