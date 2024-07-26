"""Core recipes for VASP."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.calculators.vasp import Vasp
from quacc.runners.ase import Runner
from quacc.runners.thermo import ThermoRunner
from quacc.schemas.ase import summarize_vib_and_thermo
from quacc.schemas.vasp import summarize_vasp_opt_run, vasp_summarize_run
from quacc.utils.dicts import recursive_dict_merge

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms

    from quacc.types import (
        Filenames,
        OptParams,
        SourceDirectory,
        VaspASEOptSchema,
        VaspSchema,
        VibThermoSchema,
    )


def run_and_summarize(
    atoms: Atoms,
    preset: str | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    report_mp_corrections: bool = False,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> VaspSchema:
    """
    Base job function for VASP recipes.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    calc_defaults
        Default parameters for the recipe.
    calc_swaps
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].
    report_mp_corrections
        Whether to report the Materials Project corrections in the results.
    additional_fields
        Additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    VaspSchema
        Dictionary of results
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    calc = Vasp(atoms, preset=preset, **calc_flags)
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc()

    return vasp_summarize_run(
        final_atoms,
        report_mp_corrections=report_mp_corrections,
        additional_fields=additional_fields,
    )


def run_and_summarize_opt(
    atoms: Atoms,
    preset: str | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    opt_defaults: dict[str, Any] | None = None,
    opt_params: OptParams | None = None,
    report_mp_corrections: bool = False,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> VaspASEOptSchema:
    """
    Base job function for VASP recipes with ASE optimizers.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    calc_defaults
        Default parameters for the recipe.
    calc_swaps
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].
    opt_defaults
        Default arguments for the ASE optimizer.
    opt_params
        Dictionary of custom kwargs for [quacc.runners.ase.Runner.run_opt][]
    report_mp_corrections
        Whether to report the Materials Project corrections in the results.
    additional_fields
        Additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    VaspASEOptSchema
        Dictionary of results
    """
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    calc = Vasp(atoms, preset=preset, **calc_flags)
    dyn = Runner(atoms, calc, copy_files=copy_files).run_opt(**opt_flags)

    return summarize_vasp_opt_run(
        dyn,
        report_mp_corrections=report_mp_corrections,
        additional_fields=additional_fields,
    )


def run_and_summarize_vib_and_thermo(
    atoms: Atoms,
    energy: float = 0.0,
    temperature: float = 298.15,
    pressure: float = 1.0,
    thermo_method: Literal["ideal_gas", "harmonic"] = "ideal_gas",
    preset: str | None = None,
    calc_defaults: dict[str, Any] | None = None,
    calc_swaps: dict[str, Any] | None = None,
    vib_kwargs: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> VibThermoSchema:
    """
    Base job function for VASP recipes with ASE vibrational analysis.

    Parameters
    ----------
    atoms
        Atoms object
    energy
        Energy of the system
    temperature
        Temperature of the system
    pressure
        Pressure of the system
    thermo_method
        Method to use for thermochemistry. Options are "harmonic" or "ideal_gas".
    preset
        Preset to use from `quacc.calculators.vasp.presets`.
    calc_defaults
        Default parameters for the recipe.
    calc_swaps
        Dictionary of custom kwargs for the Vasp calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.calculators.vasp.vasp.Vasp][].
    vib_kwargs
        Dictionary of custom kwargs for [quacc.runners.ase.Runner.run_vib][]
    additional_fields
        Additional fields to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    VibThermoSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_vib_and_thermo][].
        See the type-hint for the data structure.
    """

    # Set defaults
    calc_flags = recursive_dict_merge(calc_defaults, calc_swaps)

    calc = Vasp(atoms, preset=preset, **calc_flags)
    vibrations = Runner(atoms, calc, copy_files=copy_files).run_vib(
        vib_kwargs=vib_kwargs
    )
    thermo_runner = ThermoRunner(atoms, vibrations.get_frequencies(), energy=energy)
    if thermo_method == "harmonic":
        thermo_object = thermo_runner.run_harmonic_thermo()
    elif thermo_method == "ideal_gas":
        thermo_object = thermo_runner.run_ideal_gas()

    return summarize_vib_and_thermo(
        vibrations,
        thermo_object,
        atoms=atoms,
        temperature=temperature,
        pressure=pressure,
        additional_fields=additional_fields,
    )
