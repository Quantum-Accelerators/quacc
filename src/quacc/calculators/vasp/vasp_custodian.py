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

if TYPE_CHECKING:
    from typing import Callable, TypedDict

    class VaspJobKwargs(TypedDict, total=False):
        """
        Type hint for `vasp_job_kwargs` in in [quacc.calculators.vasp.vasp_custodian.run_custodian][].
        """

        output_file: str  # default = "vasp.out"
        stderr_file: str  # default = "std_err.txt"
        suffix: str  # default = ""
        final: bool  # default = True
        backup: bool  # default = True
        auto_npar: bool  # default = False
        auto_gamma: bool  # default = True
        settings_override: dict | None  # default = None
        copy_magmom: bool  # default = False
        auto_continue: bool  # default = False

    class CustodianKwargs(TypedDict, total=False):
        """
        Type hint for `custodian_kwargs` in [quacc.calculators.vasp.vasp_custodian.run_custodian][].
        """

        max_errors_per_job: int | None  # default = None
        polling_time_step: int  # default = 10
        monitor_freq: int  # default = 10
        skip_over_errors: bool  # default = False
        gzipped_output: bool  # default = False
        checkpoint: bool  # default = False
        terminate_func: Callable | None  # default = None
        terminate_on_nonzero_returncode: bool  # default = False


def run_custodian(
    vasp_parallel_cmd: str | None = None,
    vasp_cmd: str | None = None,
    vasp_gamma_cmd: str | None = None,
    vasp_custodian_max_errors: int | None = None,
    vasp_custodian_wall_time: float | None = None,
    vtst_fixes: bool | None = None,
    vasp_custodian_handlers: list[str] | None = None,
    vasp_custodian_validators: list[str] | None = None,
    scratch_dir: str | None = None,
    directory: str | None = "./",
    vasp_job_kwargs: VaspJobKwargs | None = None,
    custodian_kwargs: CustodianKwargs | None = None,
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

    from quacc import SETTINGS

    # Set defaults
    vasp_parallel_cmd = os.path.expandvars(
        SETTINGS.VASP_PARALLEL_CMD if vasp_parallel_cmd is None else vasp_parallel_cmd
    )
    vasp_cmd = SETTINGS.VASP_CMD if vasp_cmd is None else vasp_cmd
    vasp_gamma_cmd = (
        SETTINGS.VASP_GAMMA_CMD if vasp_gamma_cmd is None else vasp_gamma_cmd
    )
    vasp_custodian_max_errors = (
        SETTINGS.VASP_CUSTODIAN_MAX_ERRORS
        if vasp_custodian_max_errors is None
        else vasp_custodian_max_errors
    )
    vasp_custodian_wall_time = (
        SETTINGS.VASP_CUSTODIAN_WALL_TIME
        if vasp_custodian_wall_time is None
        else vasp_custodian_wall_time
    )
    vtst_fixes = SETTINGS.VASP_CUSTODIAN_VTST if vtst_fixes is None else vtst_fixes
    vasp_custodian_handlers = (
        SETTINGS.VASP_CUSTODIAN_HANDLERS
        if vasp_custodian_handlers is None
        else vasp_custodian_handlers
    )
    vasp_custodian_validators = (
        SETTINGS.VASP_CUSTODIAN_VALIDATORS
        if vasp_custodian_validators is None
        else vasp_custodian_validators
    )

    # Handlers for VASP
    handlers = []
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
    for handler_flag in vasp_custodian_handlers:
        if handler_flag not in handlers_dict:
            msg = f"Unknown VASP error handler: {handler_flag}"
            raise ValueError(msg)
        handlers.append(handlers_dict[handler_flag])

    validators = []
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
