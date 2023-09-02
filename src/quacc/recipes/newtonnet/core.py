"""
Core recipes for the NewtonNet code
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE
from ase.vibrations.data import VibrationsData
from monty.dev import requires

from quacc import SETTINGS, job
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import (
    summarize_opt_run,
    summarize_run,
    summarize_thermo,
    summarize_vib_run,
)
from quacc.utils.calc import run_ase_opt, run_calc
from quacc.utils.dicts import merge_dicts
from quacc.utils.thermo import ideal_gas

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

    class FreqSchema(VibSchema):
        thermo: ThermoSchema


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def static_job(
    atoms: Atoms | dict,
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
    calc_swaps
        Dictionary of custom kwargs for the newtonnet calculator. Overrides the
        following defaults:

        ```python
        {
            "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
            "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        }
        ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results, specified in `quacc.schemas.ase.RunSchema`
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

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
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def relax_job(
    atoms: Atoms | dict,
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
    calc_swaps
        Dictionary of custom kwargs for the newtonnet calculator. Overrides the
        following defaults:

        ```python
        {
            "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
            "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        }
        ```
    opt_swaps
        Optional swaps for the optimization parameters. Overrides the following
        defaults:

        ```python
        {"fmax": 0.01, "max_steps": 1000, "optimizer": Sella or FIRE}
        ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    OptSchema
        Dictionary of results, specified in `quacc.schemas.ase.OptSchema``
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}
    opt_swaps = opt_swaps or {}

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": Sella or FIRE}

    flags = merge_dicts(defaults, calc_swaps)
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    if "sella.optimize" in opt_flags.get("optimizer", FIRE).__module__:
        opt_flags["order"] = 0

    atoms.calc = NewtonNet(**flags)
    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_swaps)

    return _add_stdev_and_hess(
        summarize_opt_run(dyn, additional_fields={"name": "NewtonNet Relax"})
    )


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def freq_job(
    atoms: Atoms | dict,
    temperature: float = 298.15,
    pressure: float = 1.0,
    calc_swaps: dict | None = None,
) -> FreqSchema:
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
    calc_swaps
        Optional swaps for the calculator. Overrides the following
        defaults:

        ```python
        {
            "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
            "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        }
        ```

    Returns
    -------
    FreqSchema
        Dictionary of results
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    flags = merge_dicts(defaults, calc_swaps)

    ml_calculator = NewtonNet(**flags)
    atoms.calc = ml_calculator

    ml_calculator.calculate(atoms)
    hessian = ml_calculator.results["hessian"]
    vib = VibrationsData(atoms, hessian)
    vib_summary = summarize_vib_run(
        vib, additional_fields={"name": "NewtonNet Frequency Analysis"}
    )

    igt = ideal_gas(
        atoms, vib.get_frequencies(), energy=ml_calculator.results["energy"]
    )
    vib_summary["thermo"] = summarize_thermo(
        igt, temperature=temperature, pressure=pressure
    )

    return vib_summary


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
        ml_calculator.calculate(conf["atoms"])
        conf["hessian"] = ml_calculator.results["hessian"]
        conf["energy_std"] = ml_calculator.results["energy_disagreement"]
        conf["forces_std"] = ml_calculator.results["forces_disagreement"]
        conf["hessian_std"] = ml_calculator.results["hessian_disagreement"]

    return summary
