"""Core recipes for the NewtonNet code."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE
from ase.vibrations.data import VibrationsData
from monty.dev import requires

from quacc import SETTINGS, job
from quacc.builders.thermo import build_ideal_gas
from quacc.runners.ase import run_calc, run_opt
from quacc.schemas.ase import (
    summarize_ideal_gas_thermo,
    summarize_opt_run,
    summarize_run,
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
    from typing import Any

    from ase import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema, ThermoSchema, VibSchema

    class FreqSchema(RunSchema):
        vib: VibSchema
        thermo: ThermoSchema


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
def static_job(
    atoms: Atoms, copy_files: list[str] | None = None, **kwargs
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    copy_files
        Files to copy to the runtime directory.
    **kwargs
        Custom kwargs for the NewtonNet calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

        !!! Info "Calculator defaults"

            ```python
            {
                "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
                "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
            }
            ```

    Returns
    -------
    RunSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_run][]
    """

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    flags = merge_dicts(defaults, kwargs)

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
    atoms: Atoms,
    opt_params: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
    **kwargs,
) -> OptSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].

        !!! Info "Optimizer defaults"

            ```python
            {"fmax": 0.01, "max_steps": 1000, "optimizer": Sella or FIRE}
            ```
    copy_files
        Files to copy to the runtime directory.
    **kwargs
        Dictionary of custom kwargs for the NewtonNet calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

        !!! Info "Calculator defaults"

            ```python
            {
                "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
                "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
            }
            ```

    Returns
    -------
    OptSchema
        Dictionary of results, specified in [quacc.schemas.ase.summarize_opt_run][]
    """

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": Sella or FIRE}

    flags = merge_dicts(defaults, kwargs)
    opt_flags = merge_dicts(opt_defaults, opt_params)

    atoms.calc = NewtonNet(**flags)
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
    copy_files: list[str] | None = None,
    **kwargs,
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
    copy_files
        Files to copy to the runtime directory.
    **kwargs
        Custom kwargs for the NewtonNet calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

        !!! Info "Calculator defaults"

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

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    flags = merge_dicts(defaults, kwargs)

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

    igt = build_ideal_gas(final_atoms, vib.get_frequencies(), energy=energy)
    summary["thermo"] = summarize_ideal_gas_thermo(
        igt,
        temperature=temperature,
        pressure=pressure,
        additional_fields={"name": "ASE Thermo Analysis"},
    )

    return summary


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
