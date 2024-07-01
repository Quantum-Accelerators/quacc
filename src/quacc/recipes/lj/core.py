"""
Core recipes for Lennard-Jones Potential.

NOTE: This set of minimal recipes is mainly for demonstration purposes
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.lj import LennardJones

from quacc import job
from quacc.runners.ase import Runner
from quacc.runners.thermo import ThermoRunner
from quacc.schemas.ase import summarize_opt_run, summarize_run, summarize_vib_and_thermo

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.runners.ase import OptParams, VibKwargs
    from quacc.schemas._aliases.ase import OptSchema, RunSchema, VibThermoSchema
    from quacc.types import Filenames, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a static calculation.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Dictionary of custom kwargs for the LJ calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.lj.LennardJones] calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """
    calc = LennardJones(**calc_kwargs)
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    return summarize_run(final_atoms, atoms, additional_fields={"name": "LJ Static"})


@job
def relax_job(
    atoms: Atoms,
    opt_params: OptParams | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the LJ calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.lj.LennardJones] calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """
    opt_params = opt_params or {}

    calc = LennardJones(**calc_kwargs)
    dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(**opt_params)

    return summarize_opt_run(dyn, additional_fields={"name": "LJ Relax"})


@job
def freq_job(
    atoms: Atoms,
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    vib_kwargs: VibKwargs | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Dictionary of custom kwargs for the LJ calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the [ase.calculators.lj.LennardJones] calculator.

    Returns
    -------
    VibThermoSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_vib_and_thermo][].
        See the type-hint for the data structure.
    """
    vib_kwargs = vib_kwargs or {}

    calc = LennardJones(**calc_kwargs)
    vibrations = Runner(atoms, calc, copy_files=copy_files).run_vib(
        vib_kwargs=vib_kwargs
    )
    igt = ThermoRunner(
        atoms, vibrations.get_frequencies(), energy=energy
    ).run_ideal_gas()

    return summarize_vib_and_thermo(
        vibrations,
        igt,
        temperature=temperature,
        pressure=pressure,
        additional_fields={"name": "LJ Frequency and Thermo"},
    )
