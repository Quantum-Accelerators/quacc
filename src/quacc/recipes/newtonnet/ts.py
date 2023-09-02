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
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import summarize_opt_run
from quacc.utils.calc import run_ase_opt
from quacc.utils.dicts import merge_dicts

try:
    from sella import IRC, Sella
except ImportError:
    Sella = None

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None

if TYPE_CHECKING:
    from typing import Literal

    import numpy as np
    from ase import Atoms

    from quacc.recipes.newtonnet.core import FreqSchema
    from quacc.schemas.ase import OptSchema
    from quacc.utils.wflows import Job

    class TSSchema(OptSchema):
        freq: FreqSchema | None

    class IRCSchema(OptSchema):
        freq: FreqSchema | None

    class QuasiIRCSchema(OptSchema):
        irc: IRCSchema
        freq: FreqSchema | None


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
@requires(Sella, "Sella must be installed. Try pip install quacc[optimizers]")
def ts_job(
    atoms: Atoms | dict,
    use_custom_hessian: bool = False,
    freq_job: Job | None = _freq_job,
    freq_job_kwargs: dict | None = None,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
) -> TSSchema:
    """
    Perform a transition state (TS) job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    use_custom_hessian
        Whether to use a custom Hessian matrix.
    freq_job
        Default Job to use for the frequency analysis.
    freq_job_kwargs
        Keyword arguments to use for the `freq_job`.
    calc_swaps
        Optional swaps for the NewtonNet calculator. Overrides the
        following defaults:

        ```python
        {
            "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
            "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        }
        ```
    opt_swaps
        Optional swaps for the optimization parameters. Overrides the
        following defaults:

        ```python
        {
            "fmax": 0.01,
            "max_steps": 1000,
            "optimizer": Sella,
            "optimizer_kwargs": {"diag_every_n": 0} if use_custom_hessian else {},
        }
        ```

    Returns
    -------
    TSSchema
        Dictionary of results, specified in `quacc.schemas.ase.RunSchema`
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
        summarize_opt_run(dyn, additional_fields={"name": "NewtonNet TS"})
    )

    # Run a frequency calculation
    freq_summary = (
        freq_job.__wrapped__(opt_ts_summary, **freq_job_kwargs) if freq_job else None
    )
    opt_ts_summary["freq"] = freq_summary

    return opt_ts_summary


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
@requires(Sella, "Sella must be installed. Try pip install quacc[optimizers]")
def irc_job(
    atoms: Atoms | dict,
    direction: Literal["forward", "reverse"] = "forward",
    freq_job: Job | None = _freq_job,
    freq_job_kwargs: dict | None = None,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
) -> IRCSchema:
    """
    Perform an intrinsic reaction coordinate (IRC) job using the given atoms
    object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    direction
        The direction of the IRC calculation ("forward" or "reverse").
    freq_job
        Default Job to use for the frequency analysis.
    freq_job_kwargs
        Keyword arguments for the `freq_job`.
    calc_swaps
        Optional swaps for the calculator. Overrides the following
        defaults:

        ```python
        {
            "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
            "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        }
        ```
    opt_swaps
        Optional swaps for the optimization parameters. Overrides the
        following defaults:

        ```python
        {
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
        ```

    Returns
    -------
    IRCSchema
        A dictionary containing the IRC summary and thermodynamic summary.
    """
    atoms = fetch_atoms(atoms)
    calc_swaps = calc_swaps or {}
    opt_swaps = opt_swaps or {}
    freq_job_kwargs = freq_job_kwargs or {}
    default_settings = SETTINGS.copy()

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
    SETTINGS.CHECK_CONVERGENCE = False
    dyn = run_ase_opt(atoms, **opt_flags)
    opt_irc_summary = _add_stdev_and_hess(
        summarize_opt_run(
            dyn, additional_fields={"name": f"NewtonNet IRC: {direction}"}
        )
    )
    SETTINGS.CHECK_CONVERGENCE = default_settings.CHECK_CONVERGENCE

    # Run frequency job
    freq_summary = (
        freq_job.__wrapped__(opt_irc_summary, **freq_job_kwargs) if freq_job else None
    )
    opt_irc_summary["freq"] = freq_summary

    return opt_irc_summary


@job
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
@requires(Sella, "Sella must be installed. Try pip install quacc[optimizers]")
def quasi_irc_job(
    atoms: Atoms | dict,
    direction: Literal["forward", "reverse"] = "forward",
    freq_job: Job | None = _freq_job,
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
    freq_job
        Default Job to use for the frequency analysis.
    freq_job_kwargs
        Keyword arguments for `freq_job`.
    irc_swaps
        Optional swaps for the IRC optimization parameters. Overrides
        the following defaults: `{"max_steps": 5}`.
    opt_swaps
        Optional swaps for the optimization parameters. Overrides
        the following defaults: `{}`.

    Returns
    -------
    QuasiIRCSchema
        A dictionary containing the IRC summary, optimization summary, and
        thermodynamic summary.
    """
    irc_swaps = irc_swaps or {}
    opt_swaps = opt_swaps or {}
    freq_job_kwargs = freq_job_kwargs or {}

    irc_defaults = {"max_steps": 5}
    irc_flags = merge_dicts(irc_defaults, irc_swaps)

    # Run IRC
    irc_summary = irc_job.__wrapped__(
        atoms, direction=direction, opt_swaps=irc_flags, freq_job=None
    )

    # Run opt
    relax_summary = relax_job.__wrapped__(irc_summary, **opt_swaps)

    # Run frequency
    freq_summary = (
        freq_job.__wrapped__(relax_summary, **freq_job_kwargs) if freq_job else None
    )
    relax_summary["freq"] = freq_summary
    relax_summary["irc"] = irc_summary

    return relax_summary


def _get_hessian(atoms: Atoms) -> np.ndarray:
    """
    Calculate and retrieve the Hessian matrix for the given molecular
    configuration.

    This function takes an ASE Atoms object representing a molecular
    configuration and uses the NewtonNet machine learning calculator to
    calculate the Hessian matrix. The calculated Hessian matrix is then
    returned.

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
