"""
Core recipes for the NewtonNet code
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE
from ase.vibrations.data import VibrationsData
from monty.dev import requires

from quacc import SETTINGS, job
from quacc.runners.calc import run_ase_opt, run_calc
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
    from sella import Sella
except ImportError:
    Sella = None

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import OptSchema, RunSchema, ThermoSchema, VibSchema

    class FreqSchema(RunSchema):
        vib: VibSchema
        thermo: ThermoSchema


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
def static_job(
    atoms: Atoms | dict,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    ??? Note

        Calculator Defaults:

        ```python
        {
            "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
            "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        }
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    calc_swaps
        Dictionary of custom kwargs for the newtonnet calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """
    atoms = fetch_atoms(atoms)

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    flags = merge_dicts(defaults, calc_swaps)

    atoms.calc = NewtonNet(**flags)
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "NewtonNet Static"},
    )


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
def relax_job(
    atoms: Atoms | dict,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> OptSchema:
    """
    Relax a structure.

    ??? Note

        Calculator Defaults:

        ```python
        {
            "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
            "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        }
        ```

        Optimizer Defaults:

        ```python
        {"fmax": 0.01, "max_steps": 1000, "optimizer": Sella or FIRE}
        ```

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    calc_swaps
        Dictionary of custom kwargs for the newtonnet calculator.
    opt_swaps
        Optional swaps for the optimization parameters.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_opt_run][]
    """
    atoms = fetch_atoms(atoms)

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": Sella or FIRE}

    flags = merge_dicts(defaults, calc_swaps)
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    atoms.calc = NewtonNet(**flags)
    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_flags)

    return _add_stdev_and_hess(
        summarize_opt_run(dyn, additional_fields={"name": "NewtonNet Relax"})
    )


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
def freq_job(
    atoms: Atoms | dict,
    temperature: float = 298.15,
    pressure: float = 1.0,
    calc_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> FreqSchema:
    """
    Perform a frequency calculation using the given atoms object.

    ??? Note

        Calculator Defaults:

        ```python
        {
            "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
            "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        }
        ```

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    temperature
        The temperature for the thermodynamic analysis.
    pressure
        The pressure for the thermodynamic analysis.
    calc_swaps
        Optional swaps for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    FreqSchema
        Dictionary of results
    """
    atoms = fetch_atoms(atoms)

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    flags = merge_dicts(defaults, calc_swaps)

    ml_calculator = NewtonNet(**flags)
    atoms.calc = ml_calculator
    final_atoms = run_calc(atoms, copy_files=copy_files)

    summary = summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "NewtonNet Hessian"},
    )
    energy = summary["results"]["energy"]
    hessian = summary["results"]["hessian"]

    vib = VibrationsData(final_atoms, hessian)
    summary["vib"] = summarize_vib_run(
        vib, additional_fields={"name": "ASE Vibrations Analysis"}
    )

    igt = ideal_gas(final_atoms, vib.get_frequencies(), energy=energy)
    summary["thermo"] = summarize_thermo(
        igt,
        temperature=temperature,
        pressure=pressure,
        additional_fields={"name": "ASE Thermo Analysis"},
    )

    return summary


def _add_stdev_and_hess(summary: dict) -> dict:
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

    for conf in summary["trajectory"]:
        ml_calculator = NewtonNet(
            model_path=SETTINGS.NEWTONNET_MODEL_PATH,
            settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        )
        atoms = conf["atoms"]
        atoms.calc = ml_calculator
        results = run_calc(atoms).calc.results
        conf["hessian"] = results["hessian"]
        conf["energy_std"] = results["energy_disagreement"]
        conf["forces_std"] = results["forces_disagreement"]
        conf["hessian_std"] = results["hessian_disagreement"]

    return summary
