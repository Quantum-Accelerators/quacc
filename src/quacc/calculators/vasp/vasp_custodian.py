"""Custodian handlers for VASP."""

from __future__ import annotations

import os
import shlex
from typing import TYPE_CHECKING

from custodian import Custodian
from custodian.vasp.handlers import (
    FrozenJobErrorHandler,
    IncorrectSmearingHandler,
    KspacingMetalHandler,
    LargeSigmaHandler,
    MeshSymmetryErrorHandler,
    NonConvergingErrorHandler,
    PositiveEnergyErrorHandler,
    PotimErrorHandler,
    StdErrHandler,
    UnconvergedErrorHandler,
    VaspErrorHandler,
    WalltimeHandler,
)
from custodian.vasp.jobs import VaspJob
from custodian.vasp.validators import VaspFilesValidator, VasprunXMLValidator

from quacc import QuaccDefault, get_settings

if TYPE_CHECKING:
    from pathlib import Path

    from src.quacc.settings import QuaccSettings

    from quacc.types import DefaultSetting, VaspCustodianKwargs, VaspJobKwargs


def run_custodian(
    vasp_parallel_cmd: str | DefaultSetting = QuaccDefault,
    vasp_cmd: str | DefaultSetting = QuaccDefault,
    vasp_gamma_cmd: str | DefaultSetting = QuaccDefault,
    vasp_custodian_max_errors: int | DefaultSetting = QuaccDefault,
    vasp_custodian_wall_time: float | DefaultSetting = QuaccDefault,
    vtst_fixes: bool | DefaultSetting = QuaccDefault,
    vasp_custodian_handlers: list[str] | None | DefaultSetting = QuaccDefault,
    vasp_custodian_validators: list[str] | None | DefaultSetting = QuaccDefault,
    scratch_dir: str | None = None,
    directory: str | Path | None = None,
    vasp_job_kwargs: VaspJobKwargs | None = None,
    custodian_kwargs: VaspCustodianKwargs | None = None,
) -> list[list[dict]]:
    """
    Function to run VASP Custodian.

    Parameters
    ----------
    vasp_parallel_cmd
        VASP parallel command, e.g. "srun -N 2 --ntasks-per-node=24". Defaults
        to the $VASP_PARALLEL_CMD environment variable in settings.
    vasp_cmd
        VASP command. Defaults to "vasp_std" in settings.
    vasp_gamma_cmd
        VASP gamma command. Defaults to "vasp_gam" in settings.
    vasp_custodian_max_errors
        Maximum number of errors to allow before stopping the run. Defaults to 5
        in settings.
    vasp_custodian_wall_time
        Maximum wall time to allow before creating a STOPCAR. Defaults to
        infinity in settings.
    vtst_fixes
        Whether to apply VTST input swaps. Defaults to False in settings.
    vasp_custodian_handlers
        List of handlers to use in Custodian. See settings for list.
    vasp_custodian_validators
        List of validators to use in Custodian. See settings for list.
    scratch_dir
        Scratch directory to use. Defaults to None.
    directory
        Directory to run the calculation in. Defaults to None.
    vasp_job_kwargs
        Keyword arguments to pass to the Custodian VaspJob. Defaults to None.
    custodian_kwargs
        Any remaining keyword arguments to pass to Custodian. Defaults to None.

    Returns
    -------
    list[list[dict]]
        List of errors from each Custodian job.
    """
    # Adapted from atomate2.vasp.run.run_vasp
    settings: QuaccSettings = get_settings()

    # Set defaults
    vasp_parallel_cmd = os.path.expandvars(
        settings.VASP_PARALLEL_CMD
        if vasp_parallel_cmd == QuaccDefault
        else vasp_parallel_cmd
    )
    vasp_cmd = settings.VASP_CMD if vasp_cmd == QuaccDefault else vasp_cmd
    vasp_gamma_cmd = (
        settings.VASP_GAMMA_CMD if vasp_gamma_cmd == QuaccDefault else vasp_gamma_cmd
    )
    vasp_custodian_max_errors = (
        settings.VASP_CUSTODIAN_MAX_ERRORS
        if vasp_custodian_max_errors == QuaccDefault
        else vasp_custodian_max_errors
    )
    vasp_custodian_wall_time = (
        settings.VASP_CUSTODIAN_WALL_TIME
        if vasp_custodian_wall_time == QuaccDefault
        else vasp_custodian_wall_time
    )
    vtst_fixes = (
        settings.VASP_CUSTODIAN_VTST if vtst_fixes == QuaccDefault else vtst_fixes
    )
    vasp_custodian_handlers = (
        settings.VASP_CUSTODIAN_HANDLERS
        if vasp_custodian_handlers == QuaccDefault
        else vasp_custodian_handlers
    )

    vasp_custodian_validators = (
        settings.VASP_CUSTODIAN_VALIDATORS
        if vasp_custodian_validators == QuaccDefault
        else vasp_custodian_validators
    )

    # Handlers for VASP
    handlers_dict = {
        "VaspErrorHandler": VaspErrorHandler(vtst_fixes=vtst_fixes),
        "FrozenJobErrorHandler": FrozenJobErrorHandler(),
        "IncorrectSmearingHandler": IncorrectSmearingHandler(),
        "LargeSigmaHandler": LargeSigmaHandler(),
        "MeshSymmetryErrorHandler": MeshSymmetryErrorHandler(),
        "NonConvergingErrorHandler": NonConvergingErrorHandler(),
        "PositiveEnergyErrorHandler": PositiveEnergyErrorHandler(),
        "PotimErrorHandler": PotimErrorHandler(),
        "StdErrHandler": StdErrHandler(),
        "UnconvergedErrorHandler": UnconvergedErrorHandler(),
        "WalltimeHandler": WalltimeHandler(),
        "KspacingMetalHandler": KspacingMetalHandler(),
    }
    validators_dict = {
        "VaspFilesValidator": VaspFilesValidator(),
        "VasprunXMLValidator": VasprunXMLValidator(),
    }

    handlers = []
    if vasp_custodian_handlers is None:
        vasp_custodian_handlers = []

    for handler_flag in vasp_custodian_handlers:
        if handler_flag not in handlers_dict:
            msg = f"Unknown VASP error handler: {handler_flag}"
            raise ValueError(msg)
        handlers.append(handlers_dict[handler_flag])

    validators = []
    if vasp_custodian_validators is None:
        vasp_custodian_validators = []
    for validator_flag in vasp_custodian_validators:
        if validator_flag not in validators_dict:
            msg = f"Unknown VASP validator: {validator_flag}"
            raise ValueError(msg)
        validators.append(validators_dict[validator_flag])

    # Populate settings
    full_vasp_cmd = f"{vasp_parallel_cmd} {vasp_cmd}"
    full_vasp_gamma_cmd = f"{vasp_parallel_cmd} {vasp_gamma_cmd}"

    # Run VASP
    vasp_job_kwargs = {} if vasp_job_kwargs is None else vasp_job_kwargs
    custodian_kwargs = {} if custodian_kwargs is None else custodian_kwargs
    split_vasp_cmd = shlex.split(full_vasp_cmd)
    split_vasp_gamma_cmd = shlex.split(full_vasp_gamma_cmd)
    vasp_job_kwargs["gamma_vasp_cmd"] = split_vasp_gamma_cmd

    # Run with Custodian
    jobs = [VaspJob(split_vasp_cmd, **vasp_job_kwargs)]

    if vasp_custodian_wall_time:
        handlers = [
            *list(handlers),
            WalltimeHandler(wall_time=vasp_custodian_wall_time),
        ]

    c = Custodian(
        handlers,
        jobs,
        validators=validators,
        max_errors=vasp_custodian_max_errors,
        scratch_dir=scratch_dir,
        directory=directory,
        **custodian_kwargs,
    )

    return c.run()


if __name__ == "__main__":
    run_custodian()
