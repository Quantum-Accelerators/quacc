"""Core recipes for the NewtonNet code."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.vibrations.data import VibrationsData
from monty.dev import requires

from quacc import SETTINGS, job
from quacc.runners.ase import run_calc, run_opt
from quacc.runners.thermo import run_ideal_gas
from quacc.schemas.ase import (
    summarize_opt_run,
    summarize_run,
    summarize_vib_and_thermo,
)
from quacc.utils.dicts import recursive_dict_merge

try:
    from sella import Sella
except ImportError:
    Sella = None

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema, VibThermoSchema
    from quacc.utils.files import Filenames, SourceDirectory


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
def static_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)

    atoms.calc = NewtonNet(**calc_flags)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms, atoms, additional_fields={"name": "NewtonNet Static"}
    )


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
def relax_job(
    atoms: Atoms,
    opt_params: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> OptSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Dictionary of custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_opt_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {"optimizer": Sella} if Sella else {}

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    atoms.calc = NewtonNet(**calc_flags)
    dyn = run_opt(atoms, copy_files=copy_files, **opt_flags)

    return _add_stdev_and_hess(
        summarize_opt_run(dyn, additional_fields={"name": "NewtonNet Relax"})
    )


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
def freq_job(
    atoms: Atoms,
    temperature: float = 298.15,
    pressure: float = 1.0,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> VibThermoSchema:
    """
    Perform a frequency calculation using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    temperature
        The temperature for the thermodynamic analysis.
    pressure
        The pressure for the thermodynamic analysis.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    VibThermoSchema
        Dictionary of results. See the type-hint for the data structure.
    """

    calc_defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        "hess_method": "autograd",
    }
    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)

    ml_calculator = NewtonNet(**calc_flags)
    atoms.calc = ml_calculator
    final_atoms = run_calc(atoms, copy_files=copy_files)

    summary = summarize_run(
        final_atoms, atoms, additional_fields={"name": "NewtonNet Hessian"}
    )
    energy = summary["results"]["energy"]
    hessian = summary["results"]["hessian"]

    vib = VibrationsData(final_atoms, hessian)

    igt = run_ideal_gas(final_atoms, vib.get_frequencies(), energy=energy)

    return summarize_vib_and_thermo(
        vib,
        igt,
        temperature=temperature,
        pressure=pressure,
        additional_fields={"name": "ASE Vibrations and Thermo Analysis"},
    )


def _add_stdev_and_hess(summary: dict[str, Any]) -> dict[str, Any]:
    """
    Calculate and add standard deviation values and Hessians to the summary.

    This function takes a summary dictionary containing information about a
    molecular trajectory and calculates the standard deviation of various
    properties using the NewtonNet machine learning calculator. It adds the
    calculated standard deviation values and Hessians to each configuration in
    the trajectory.

    Parameters
    ----------
    summary
        A dictionary containing information about the molecular trajectory.

    Returns
    -------
    Dict
        The modified summary dictionary with added standard deviation and
        Hessian values.
    """

    for i, atoms in enumerate(summary["trajectory"]):
        ml_calculator = NewtonNet(
            model_path=SETTINGS.NEWTONNET_MODEL_PATH,
            settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        )
        atoms.calc = ml_calculator
        results = run_calc(atoms).calc.results
        summary["trajectory_results"][i]["hessian"] = results["hessian"]
        summary["trajectory_results"][i]["energy_std"] = results["energy_disagreement"]
        summary["trajectory_results"][i]["forces_std"] = results["forces_disagreement"]
        summary["trajectory_results"][i]["hessian_std"] = results[
            "hessian_disagreement"
        ]

    return summary
