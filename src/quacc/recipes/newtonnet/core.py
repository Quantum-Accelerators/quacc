"""
Core recipes for the NewtonNet code
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from ase.optimize import FIRE
from ase.vibrations.data import VibrationsData
from monty.dev import requires

from quacc import SETTINGS, job
from quacc.schemas.ase import (
    summarize_opt_run,
    summarize_run,
    summarize_thermo_run,
    summarize_vib_run,
)
from quacc.utils.calc import run_ase_opt, run_calc
from quacc.utils.dicts import merge_dicts
from quacc.utils.thermo import ideal_gas
from quacc.utils.wflows import fetch_atoms

if TYPE_CHECKING:
    import numpy as np
    from ase import Atoms
    from ase.optimize.optimize import Optimizer

    from quacc.schemas.ase import OptSchema, RunSchema, ThermoSchema, VibSchema


try:
    from sella import IRC, Sella
except ImportError:
    Sella = None

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None


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
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    calc_swaps
        Dictionary of custom kwargs for the newtonnet calculator
    opt_swaps
        Optional swaps for the optimization parameters
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        A dictionary containing the results of the calculation.
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}

    # Define calculator
    atoms.calc = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        **calc_swaps,
    )
    final_atoms = run_calc(atoms, copy_files=copy_files)

    return summarize_run(
        final_atoms, input_atoms=atoms, additional_fields={"name": "NewtonNet S"}
    )


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def relax_job(
    atoms: Atoms | dict,
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: Optimizer = Sella or FIRE,
    newtonnet_kwargs: dict | None = None,
    optimizer_kwargs: dict | None = None,
) -> OptSchema:
    """
    Relax a structure.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as the value
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        Optimizer class to use for the relaxation.
    newtonnet_kwargs
        Dictionary of custom kwargs for the newtonnet calculator.
    optimizer_kwargs
        Dictionary of kwargs for the optimizer.

    Returns
    -------
    OptSchema
        A dictionary containing the results of the calculation.
    """
    atoms = fetch_atoms(atoms)
    newtonnet_kwargs = newtonnet_kwargs or {}
    optimizer_kwargs = optimizer_kwargs or {}
    if "sella.optimize" in optimizer.__module__:
        optimizer_kwargs["order"] = 0

    atoms.calc = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
    )
    dyn = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        optimizer_kwargs=optimizer_kwargs,
    )
    return _add_stdev_and_hess(
        summarize_opt_run(dyn, additional_fields={"name": "NewtonNet Relax"})
    )


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def freq_job(
    atoms: Atoms | dict,
    temperature: float = 298.15,
    pressure: float = 1.0,
) -> dict[Literal["vib", "thermo"], VibSchema | ThermoSchema]:
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

    Returns
    -------
    dict
        Summary of the frequency calculation and thermo calculations.
    """
    atoms = fetch_atoms(atoms)

    # Define calculator
    ml_calculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
    )
    atoms.calc = ml_calculator

    # Run calculator
    ml_calculator.calculate(atoms)
    hessian = ml_calculator.results["hessian"]
    vib = VibrationsData(atoms, hessian)

    # Make IdealGasThermo object
    igt = ideal_gas(
        atoms, vib.get_frequencies(), energy=ml_calculator.results["energy"]
    )

    return {
        "vib": summarize_vib_run(
            vib, additional_fields={"name": "NewtonNet Vibrations"}
        ),
        "thermo": summarize_thermo_run(
            igt,
            temperature=temperature,
            pressure=pressure,
            additional_fields={"name": "NewtonNet Thermo"},
        ),
    }


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
@requires(Sella, "Sella must be installed. Try pip install quacc[optimizers]")
def ts_job(
    atoms: Atoms | dict,
    use_custom_hessian: bool = False,
    temperature: float = 298.15,
    pressure: float = 1.0,
    check_convergence: bool = True,
    opt_swaps: dict | None = None,
) -> dict[Literal["ts", "thermo"], OptSchema | ThermoSchema]:
    """
    Perform a transition state (TS) job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    use_custom_hessian
        Whether to use a custom Hessian matrix.
    temperature
        The temperature for the frequency calculation in Kelvins.
    pressure
        The pressure for the frequency calculation in bar.
    check_convergence
        Whether to check the convergence of the optimization.
    opt_swaps
        Optional swaps for the optimization parameters.

    Returns
    -------
    dict
        A dictionary containing the TS summary and thermodynamic summary.
    """
    atoms = fetch_atoms(atoms)
    opt_swaps = opt_swaps or {}

    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella,
        "optimizer_kwargs": {"diag_every_n": 0} if use_custom_hessian else {},
    }

    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    # Define calculator
    atoms.calc = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
    )

    if use_custom_hessian:
        if opt_flags["optimizer"].__name__ != "Sella":
            raise ValueError("Custom hessian can only be used with Sella.")

        opt_flags["optimizer_kwargs"]["hessian_function"] = _get_hessian

    ml_calculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
    )
    atoms.calc = ml_calculator

    # Run the TS optimization
    dyn = run_ase_opt(atoms, **opt_flags)

    ts_summary = summarize_opt_run(
        dyn,
        check_convergence=check_convergence,
        additional_fields={"name": "NewtonNet TS"},
    )

    ts_summary = _add_stdev_and_hess(ts_summary)

    # Run a frequency calculation
    thermo_summary = freq_job(ts_summary, temperature=temperature, pressure=pressure)

    return {"ts": ts_summary, "thermo": thermo_summary}


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
@requires(Sella, "Sella must be installed. Try pip install quacc[optimizers]")
def irc_job(
    atoms: Atoms | dict,
    fmax: float = 0.01,
    max_steps: int = 1000,
    temperature: float = 298.15,
    pressure: float = 1.0,
    check_convergence: bool = False,
    opt_swaps: dict | None = None,
) -> dict[Literal["irc", "thermo"], OptSchema | ThermoSchema]:
    """
    Perform an intrinsic reaction coordinate (IRC) job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    temperature
        The temperature for the frequency calculation in Kelvins.
    pressure
        The pressure for the frequency calculation in bar.
    check_convergence
        Whether to check the convergence of the optimization.
    opt_swaps
        Optional swaps for the optimization parameters.

    Returns
    -------
    dict
        A dictionary containing the IRC summary and thermodynamic summary.
    """
    atoms = fetch_atoms(atoms)
    opt_swaps = opt_swaps or {}

    opt_defaults = {
        "optimizer": IRC,
        "optimizer_kwargs": {
            "dx": 0.1,
            "eta": 1e-4,
            "gamma": 0.4,
            "keep_going": True,
        },
        "run_kwargs": {
            "direction": "forward",
        },
    }

    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    # Define calculator
    atoms.calc = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
    )

    # Run IRC
    dyn = run_ase_opt(atoms, fmax=fmax, max_steps=max_steps, **opt_flags)
    summary_irc = summarize_opt_run(
        dyn,
        check_convergence=check_convergence,
        additional_fields={"name": "NewtonNet IRC"},
    )

    summary_irc = _add_stdev_and_hess(summary_irc)

    # Run frequency job
    thermo_summary = freq_job.original_func(
        summary_irc, temperature=temperature, pressure=pressure
    )
    return {"irc": summary_irc, "thermo": thermo_summary}


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
@requires(Sella, "Sella must be installed. Try pip install quacc[optimizers]")
def quasi_irc_job(
    atoms: Atoms | dict,
    direction: Literal["forward", "reverse"] = "forward",
    temperature: float = 298.15,
    pressure: float = 1.0,
    irc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
) -> dict[Literal["irc", "opt", "thermo"], OptSchema | ThermoSchema]:
    """
    Perform a quasi-IRC job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    direction
        The direction of the IRC calculation ("forward" or "reverse").
    temperature
        The temperature for the frequency calculation in Kelvins.
    pressure
        The pressure for the frequency calculation in bar.
    irc_swaps
        Optional swaps for the IRC optimization parameters.
    opt_swaps
        Optional swaps for the optimization parameters.

    Returns
    -------
    dict
        A dictionary containing the IRC summary, optimization summary, and thermodynamic summary.
    """
    irc_swaps = irc_swaps or {}
    opt_swaps = opt_swaps or {}

    irc_defaults = {"run_kwargs": {"direction": direction.lower()}}
    irc_flags = merge_dicts(irc_defaults, irc_swaps)
    opt_swaps = opt_swaps or {}

    opt_defaults = {}
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    # Run IRC
    irc_summary = irc_job.original_func(atoms, max_steps=5, opt_swaps=irc_flags)

    # Run opt
    opt_summary = relax_job.original_func(irc_summary["irc"], **opt_flags)

    # Run frequency
    thermo_summary = freq_job(
        opt_summary,
        temperature=temperature,
        pressure=pressure,
    )

    return {"irc": irc_summary, "opt": opt_summary, "thermo": thermo_summary}


def _get_hessian(atoms: Atoms) -> np.ndarray:
    """
    Calculate and retrieve the Hessian matrix for the given molecular configuration.

    This function takes an ASE Atoms object representing a molecular configuration and uses the
    NewtonNet machine learning calculator to calculate the Hessian matrix. The calculated Hessian
    matrix is then returned.

    Parameters
    ----------
    atoms
        The ASE Atoms object representing the molecular configuration.

    Returns
    -------
    np.ndarray
        The calculated Hessian matrix, reshaped into a 2D array.
    """
    ml_calculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
    )
    ml_calculator.calculate(atoms)

    return ml_calculator.results["hessian"].reshape((-1, 3 * len(atoms)))


def _add_stdev_and_hess(summary: dict[str, any]) -> dict[str, any]:
    """
    Calculate and add standard deviation values and Hessians to the summary.

    This function takes a summary dictionary containing information about a molecular trajectory
    and calculates the standard deviation of various properties using the NewtonNet machine learning
    calculator. It adds the calculated standard deviation values and Hessians to each configuration
    in the trajectory.

    Parameters
    ----------
    summary
        A dictionary containing information about the molecular trajectory.

    Returns
    -------
    Dict[str, Any]
        The modified summary dictionary with added standard deviation and Hessian values.
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
