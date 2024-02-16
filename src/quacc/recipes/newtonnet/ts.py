"""Transition state recipes for the NewtonNet code."""

from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import SETTINGS, job, strip_decorator
from quacc.recipes.newtonnet.core import _add_stdev_and_hess, freq_job, relax_job
from quacc.runners.ase import run_opt
from quacc.schemas.ase import summarize_opt_run
from quacc.utils.dicts import recursive_dict_merge

try:
    from sella import IRC, Sella
except ImportError:
    Sella = None

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.atoms import Atoms
    from numpy.typing import NDArray

    from quacc.recipes.newtonnet.core import FreqSchema
    from quacc.schemas._aliases.ase import OptSchema

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
    atoms: Atoms,
    use_custom_hessian: bool = False,
    run_freq: bool = True,
    freq_job_kwargs: dict[str, Any] | None = None,
    opt_params: dict[str, Any] | None = None,
    **calc_kwargs,
) -> TSSchema:
    """
    Perform a transition state (TS) job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    use_custom_hessian
        Whether to use a custom Hessian matrix.
    run_freq
        Whether to run the frequency job.
    freq_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.freq_job][]
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    **calc_kwargs
        Dictionary of custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    TSSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    freq_job_kwargs = freq_job_kwargs or {}

    calc_defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        "hess_method": "autograd",
    }
    opt_defaults = {
        "optimizer": Sella,
        "optimizer_kwargs": (
            {"diag_every_n": 0, "order": 1} if use_custom_hessian else {"order": 1}
        ),
    }

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    atoms.calc = NewtonNet(**calc_flags)

    if use_custom_hessian:
        opt_flags["optimizer_kwargs"]["hessian_function"] = _get_hessian

    ml_calculator = NewtonNet(**calc_flags)
    atoms.calc = ml_calculator

    # Run the TS optimization
    dyn = run_opt(atoms, **opt_flags)
    opt_ts_summary = _add_stdev_and_hess(
        summarize_opt_run(dyn, additional_fields={"name": "NewtonNet TS"})
    )

    # Run a frequency calculation
    freq_summary = (
        strip_decorator(freq_job)(opt_ts_summary["atoms"], **freq_job_kwargs)
        if run_freq
        else None
    )
    opt_ts_summary["freq_job"] = freq_summary

    return opt_ts_summary


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def irc_job(
    atoms: Atoms,
    direction: Literal["forward", "reverse"] = "forward",
    run_freq: bool = True,
    freq_job_kwargs: dict[str, Any] | None = None,
    opt_params: dict[str, Any] | None = None,
    **calc_kwargs,
) -> IRCSchema:
    """
    Perform an intrinsic reaction coordinate (IRC) job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    direction
        The direction of the IRC calculation ("forward" or "reverse").
    run_freq
        Whether to run the frequency analysis.
    freq_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.freq_job][]
    opt_params
        Dictionary of custom kwargs for the optimization process. Set a value
        to `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to [quacc.runners.ase.run_opt][].
    **calc_kwargs
        Custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    IRCSchema
        A dictionary containing the IRC summary and thermodynamic summary.
        See the type-hint for the data structure.
    """
    freq_job_kwargs = freq_job_kwargs or {}
    default_settings = SETTINGS.model_copy()

    calc_defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {
        "optimizer": IRC,
        "optimizer_kwargs": {"dx": 0.1, "eta": 1e-4, "gamma": 0.4, "keep_going": True},
        "run_kwargs": {"direction": direction},
    }

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    # Define calculator
    atoms.calc = NewtonNet(**calc_flags)

    # Run IRC
    SETTINGS.CHECK_CONVERGENCE = False
    dyn = run_opt(atoms, **opt_flags)
    opt_irc_summary = _add_stdev_and_hess(
        summarize_opt_run(
            dyn, additional_fields={"name": f"NewtonNet IRC: {direction}"}
        )
    )
    SETTINGS.CHECK_CONVERGENCE = default_settings.CHECK_CONVERGENCE

    # Run frequency job
    freq_summary = (
        strip_decorator(freq_job)(opt_irc_summary["atoms"], **freq_job_kwargs)
        if run_freq
        else None
    )
    opt_irc_summary["freq_job"] = freq_summary

    return opt_irc_summary


@job
@requires(NewtonNet, "NewtonNet must be installed. Refer to the quacc documentation.")
@requires(Sella, "Sella must be installed. Refer to the quacc documentation.")
def quasi_irc_job(
    atoms: Atoms,
    direction: Literal["forward", "reverse"] = "forward",
    run_freq: bool = True,
    irc_job_kwargs: dict[str, Any] | None = None,
    relax_job_kwargs: dict[str, Any] | None = None,
    freq_job_kwargs: dict[str, Any] | None = None,
) -> QuasiIRCSchema:
    """
    Perform a quasi-IRC job using the given atoms object. The initial IRC job by default
    is run with `max_steps: 5`.

    Parameters
    ----------
    atoms
        The atoms object representing the system
    direction
        The direction of the IRC calculation
    run_freq
        Whether to run the frequency analysis
    irc_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.irc_job][]
    relax_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.core.relax_job][]
    freq_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.freq_job][]

    Returns
    -------
    QuasiIRCSchema
        A dictionary containing the IRC summary, optimization summary, and
        thermodynamic summary.
        See the type-hint for the data structure.
    """
    relax_job_kwargs = relax_job_kwargs or {}
    freq_job_kwargs = freq_job_kwargs or {}

    irc_job_defaults = {"max_steps": 5}
    irc_job_kwargs = recursive_dict_merge(irc_job_defaults, irc_job_kwargs)

    # Run IRC
    irc_summary = strip_decorator(irc_job)(
        atoms, direction=direction, run_freq=False, **irc_job_kwargs
    )

    # Run opt
    relax_summary = strip_decorator(relax_job)(irc_summary["atoms"], **relax_job_kwargs)

    # Run frequency
    freq_summary = (
        strip_decorator(freq_job)(relax_summary["atoms"], **freq_job_kwargs)
        if run_freq
        else None
    )
    relax_summary["freq_job"] = freq_summary
    relax_summary["irc_job"] = irc_summary

    return relax_summary


def _get_hessian(atoms: Atoms) -> NDArray:
    """
    Calculate and retrieve the Hessian matrix for the given molecular configuration.

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
    NDArray
        The calculated Hessian matrix, reshaped into a 2D array.
    """
    ml_calculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
    )
    ml_calculator.calculate(atoms)

    return ml_calculator.results["hessian"].reshape((-1, 3 * len(atoms)))
