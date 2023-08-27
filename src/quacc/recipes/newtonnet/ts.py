"""
Transition state recipes for the NewtonNet code
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import FIRE
from monty.dev import requires

from quacc import SETTINGS, job
from quacc.recipes.newtonnet.core import _add_stdev_and_hess
from quacc.recipes.newtonnet.core import freq_job as _freq_job
from quacc.recipes.newtonnet.core import relax_job
from quacc.schemas.ase import summarize_opt_run
from quacc.utils.calc import run_ase_opt
from quacc.utils.dicts import merge_dicts
from quacc.utils.wflows import fetch_atoms

try:
    from sella import IRC, Sella
except ImportError:
    Sella = None

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None

if TYPE_CHECKING:
    from typing import Literal, TypedDict

    import numpy as np
    from ase import Atoms

    from quacc.recipes.newtonnet.core import FreqSchema
    from quacc.schemas.ase import OptSchema

    class TSSchema(TypedDict):
        ts: OptSchema
        freq: FreqSchema | None
        atoms: Atoms

    class IRCSchema(TypedDict):
        irc: OptSchema
        freq: FreqSchema | None
        atoms: Atoms

    class QuasiIRCSchema(TypedDict):
        irc: IRCSchema
        opt: OptSchema
        freq: FreqSchema | None
        atoms: Atoms


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
@requires(Sella, "Sella must be installed. Try pip install quacc[optimizers]")
def ts_job(
    atoms: Atoms | dict,
    use_custom_hessian: bool = False,
    freq_job: callable | None = _freq_job,
    freq_job_kwargs: dict | None = None,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    check_convergence: bool = True,
) -> TSSchema:
    """
    Perform a transition state (TS) job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    use_custom_hessian
        Whether to use a custom Hessian matrix.
    calc_swaps
        Optional swaps for the calculator.
    opt_swaps
        Optional swaps for the optimization parameters.
    check_convergence
        Whether to check the convergence of the optimization.

    Returns
    -------
    dict
        A dictionary containing the TS summary and thermodynamic summary.
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}
    opt_swaps = opt_swaps or {}
    freq_job_kwargs = freq_job_kwargs or {}

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella,
        "optimizer_kwargs": {"diag_every_n": 0} if use_custom_hessian else {},
    }

    flags = merge_dicts(defaults, calc_swaps)
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    atoms.calc = NewtonNet(**flags)

    if use_custom_hessian:
        if opt_flags.get("optimizer", FIRE).__name__ != "Sella":
            raise ValueError("Custom hessian can only be used with Sella.")

        opt_flags["optimizer_kwargs"]["hessian_function"] = _get_hessian

    ml_calculator = NewtonNet(**flags)
    atoms.calc = ml_calculator

    # Run the TS optimization
    dyn = run_ase_opt(atoms, **opt_flags)
    opt_ts_summary = _add_stdev_and_hess(
        summarize_opt_run(
            dyn,
            check_convergence=check_convergence,
            additional_fields={"name": "NewtonNet TS"},
        )
    )

    # Run a frequency calculation
    freq_summary = (
        freq_job.undecorated(opt_ts_summary, **freq_job_kwargs) if freq_job else None
    )

    return {
        "atoms": opt_ts_summary["atoms"],
        "ts": opt_ts_summary,
        "freq": freq_summary,
    }


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
@requires(Sella, "Sella must be installed. Try pip install quacc[optimizers]")
def irc_job(
    atoms: Atoms | dict,
    direction: Literal["forward", "reverse"] = "forward",
    freq_job: callable | None = _freq_job,
    freq_job_kwargs: dict | None = None,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    check_convergence: bool = False,
) -> IRCSchema:
    """
    Perform an intrinsic reaction coordinate (IRC) job using the given atoms object.

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
    calc_swaps
        Optional swaps for the calculator.
    opt_swaps
        Optional swaps for the optimization parameters.
    check_convergence
        Whether to check the convergence of the optimization.

    Returns
    -------
    dict
        A dictionary containing the IRC summary and thermodynamic summary.
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}
    opt_swaps = opt_swaps or {}
    freq_job_kwargs = freq_job_kwargs or {}

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": IRC,
        "optimizer_kwargs": {
            "dx": 0.1,
            "eta": 1e-4,
            "gamma": 0.4,
            "keep_going": True,
        },
        "run_kwargs": {
            "direction": direction,
        },
    }

    flags = merge_dicts(defaults, calc_swaps)
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    # Define calculator
    atoms.calc = NewtonNet(**flags)

    # Run IRC
    dyn = run_ase_opt(atoms, **opt_flags)
    opt_irc_summary = _add_stdev_and_hess(
        summarize_opt_run(
            dyn,
            check_convergence=check_convergence,
            additional_fields={"name": f"NewtonNet IRC: {direction}"},
        )
    )

    # Run frequency job
    freq_summary = (
        freq_job.undecorated(opt_irc_summary, **freq_job_kwargs) if freq_job else None
    )
    return {
        "atoms": opt_irc_summary["atoms"],
        "irc": opt_irc_summary,
        "freq": freq_summary,
    }


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
@requires(Sella, "Sella must be installed. Try pip install quacc[optimizers]")
def quasi_irc_job(
    atoms: Atoms | dict,
    direction: Literal["forward", "reverse"] = "forward",
    freq_job: callable | None = _freq_job,
    freq_job_kwargs: dict | None = None,
    irc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
) -> QuasiIRCSchema:
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
    freq_job_kwargs = freq_job_kwargs or {}

    irc_defaults = {"max_steps": 5}
    irc_flags = merge_dicts(irc_defaults, irc_swaps)

    # Run IRC
    irc_summary = irc_job.undecorated(
        atoms, direction=direction, opt_swaps=irc_flags, freq_job=None
    )

    # Run opt
    relax_summary = relax_job.undecorated(irc_summary, **opt_swaps)

    # Run frequency
    freq_summary = (
        freq_job.undecorated(relax_summary, **freq_job_kwargs) if freq_job else None
    )

    return {
        "atoms": relax_summary["atoms"],
        "irc": irc_summary,
        "opt": relax_summary,
        "freq": freq_summary,
    }


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
