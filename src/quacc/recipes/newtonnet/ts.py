"""
Transition state recipes for the NewtonNet code
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import SETTINGS, job
from quacc.recipes.newtonnet.core import _add_stdev_and_hess, freq_job, relax_job
from quacc.runners.calc import run_ase_opt
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import summarize_opt_run
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

    from quacc.schemas.ase import FreqSchema, OptSchema

    class TSSchema(OptSchema):
        freq_job: FreqSchema | None

    class IRCSchema(OptSchema):
        freq_job: FreqSchema | None

    class QuasiIRCSchema(OptSchema):
        irc_job: IRCSchema
        freq_job: FreqSchema | None


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def ts_job(
    atoms: Atoms | dict,
    use_custom_hessian: bool = False,
    run_freq: bool = True,
    freq_job_kwargs: dict | None = None,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> TSSchema:
    """
    Perform a transition state (TS) job using the given atoms object.

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
        {
            "fmax": 0.01,
            "max_steps": 1000,
            "optimizer": Sella,
            "optimizer_kwargs": {"diag_every_n": 0, "order": 1}
            if use_custom_hessian
            else {"order": 1},
        }
        ```

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    use_custom_hessian
        Whether to use a custom Hessian matrix.
    run_freq
        Whether to run the frequency job.
    freq_job_kwargs
        Keyword arguments to use for the `freq_job`.
    calc_swaps
        Optional swaps for the NewtonNet calculator.
    opt_swaps
        Optional swaps for the optimization parameters.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    TSSchema
        Dictionary of results
    """
    atoms = fetch_atoms(atoms)
    freq_job_kwargs = freq_job_kwargs or {}

    defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": Sella,
        "optimizer_kwargs": {"diag_every_n": 0, "order": 1}
        if use_custom_hessian
        else {"order": 1},
    }

    flags = merge_dicts(defaults, calc_swaps)
    opt_flags = merge_dicts(opt_defaults, opt_swaps)

    atoms.calc = NewtonNet(**flags)

    if use_custom_hessian:
        opt_flags["optimizer_kwargs"]["hessian_function"] = _get_hessian

    ml_calculator = NewtonNet(**flags)
    atoms.calc = ml_calculator

    # Run the TS optimization
    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_flags)
    opt_ts_summary = _add_stdev_and_hess(
        summarize_opt_run(dyn, additional_fields={"name": "NewtonNet TS"})
    )

    # Run a frequency calculation
    freq_summary = (
        freq_job.__wrapped__(opt_ts_summary, **freq_job_kwargs) if run_freq else None
    )
    opt_ts_summary["freq_job"] = freq_summary

    return opt_ts_summary


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def irc_job(
    atoms: Atoms | dict,
    direction: Literal["forward", "reverse"] = "forward",
    run_freq: bool = True,
    freq_job_kwargs: dict | None = None,
    calc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> IRCSchema:
    """
    Perform an intrinsic reaction coordinate (IRC) job using the given atoms
    object.

    ??? Note

        Calculator Defaults:

        ```python
        {
            "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
            "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        }
        ```

        IRC Defaults:

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

        Optimizer Defaults:

        ```python
        {}
        ```

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    direction
        The direction of the IRC calculation ("forward" or "reverse").
    run_freq
        Whether to run the frequency analysis.
    freq_job_kwargs
        Keyword arguments for the `freq_job`.
    calc_swaps
        Optional swaps for the calculator.
    opt_swaps
        Optional swaps for the optimization parameters.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    IRCSchema
        A dictionary containing the IRC summary and thermodynamic summary.
    """
    atoms = fetch_atoms(atoms)
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
    dyn = run_ase_opt(atoms, copy_files=copy_files, **opt_flags)
    opt_irc_summary = _add_stdev_and_hess(
        summarize_opt_run(
            dyn, additional_fields={"name": f"NewtonNet IRC: {direction}"}
        )
    )
    SETTINGS.CHECK_CONVERGENCE = default_settings.CHECK_CONVERGENCE

    # Run frequency job
    freq_summary = (
        freq_job.__wrapped__(opt_irc_summary, **freq_job_kwargs) if run_freq else None
    )
    opt_irc_summary["freq_job"] = freq_summary

    return opt_irc_summary


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def quasi_irc_job(
    atoms: Atoms | dict,
    direction: Literal["forward", "reverse"] = "forward",
    run_freq: bool = True,
    irc_job_kwargs: dict | None = None,
    relax_job_kwargs: dict | None = None,
    freq_job_kwargs: dict | None = None,
    copy_files: list[str] | None = None,
) -> QuasiIRCSchema:
    """
    Perform a quasi-IRC job using the given atoms object.

    ??? Note

        IRC Defaults:

        ```python
        {"max_steps": 5}
        ```

        Optimizer Defaults:

        ```python
        {}
        ```

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    direction
        The direction of the IRC calculation ("forward" or "reverse").
    run_freq
        Whether to run the frequency analysis.
    irc_job_kwargs
        Keyword arguments for `irc_job`
    relax_job_kwargs
        Keyword arguments for `relax_job`
    freq_job_kwargs
        Keyword arguments for `freq_job`.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    QuasiIRCSchema
        A dictionary containing the IRC summary, optimization summary, and
        thermodynamic summary.
    """
    relax_job_kwargs = relax_job_kwargs or {}
    freq_job_kwargs = freq_job_kwargs or {}

    irc_job_defaults = {"calc_swaps": {"max_steps": 5}}
    irc_job_kwargs = merge_dicts(irc_job_defaults, irc_job_kwargs)

    # Run IRC
    irc_summary = irc_job.__wrapped__(
        atoms,
        direction=direction,
        run_freq=False,
        copy_files=copy_files,
        **irc_job_kwargs,
    )

    # Run opt
    relax_summary = relax_job.__wrapped__(irc_summary, **relax_job_kwargs)

    # Run frequency
    freq_summary = (
        freq_job.__wrapped__(relax_summary, **freq_job_kwargs) if run_freq else None
    )
    relax_summary["freq_job"] = freq_summary
    relax_summary["irc_job"] = irc_summary

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
