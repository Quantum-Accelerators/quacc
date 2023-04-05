"""
Custodian runner for VASP
"""
from __future__ import annotations

import os
import shlex
from typing import Dict, List

from custodian import Custodian
from custodian.vasp.handlers import (
    FrozenJobErrorHandler,
    IncorrectSmearingHandler,
    LargeSigmaHandler,
    MeshSymmetryErrorHandler,
    NonConvergingErrorHandler,
    PositiveEnergyErrorHandler,
    PotimErrorHandler,
    ScanMetalHandler,
    StdErrHandler,
    UnconvergedErrorHandler,
    VaspErrorHandler,
    WalltimeHandler,
)
from custodian.vasp.jobs import VaspJob
from custodian.vasp.validators import VaspFilesValidator, VasprunXMLValidator

from quacc import SETTINGS


def run_custodian(
    vasp_parallel_cmd: str = SETTINGS.VASP_PARALLEL_CMD,
    vasp_cmd: str = SETTINGS.VASP_CMD,
    vasp_gamma_cmd: str = SETTINGS.VASP_GAMMA_CMD,
    vasp_custodian_max_errors: int = SETTINGS.VASP_CUSTODIAN_MAX_ERRORS,
    vasp_custodian_wall_time: int = SETTINGS.VASP_CUSTODIAN_WALL_TIME,
    vtst_fixes: bool = SETTINGS.VASP_CUSTODIAN_VTST,
    vasp_custodian_handlers: List[str] = SETTINGS.VASP_CUSTODIAN_HANDLERS,
    vasp_custodian_validators: List[str] = SETTINGS.VASP_CUSTODIAN_VALIDATORS,
    scratch_dir: str = None,
    vasp_job_kwargs: Dict = None,
    custodian_kwargs: Dict = None,
):
    """
    Function to run VASP Custodian

    Parameters
    ----------
    vasp_parallel_cmd
        VASP parallel command, e.g. srun -N 2 --ntasks-per-node=24. Defaults to the $VASP_PARAlLEL_CMD
        environment variable in settings.
    vasp_cmd
        VASP command. Defaults to vasp_std in settings.
    vasp_gamma_cmd
        VASP gamma command. Defaults to vasp_gam in settings.
    vasp_custodian_max_errors
        Maximum number of errors to allow before stopping the run. Defaults to 5 in settings.
    vasp_custodian_wall_time
        Maximum wall time to allow before creating a STOPCAR. Defaults to None in settings.
    vtst_fixes
        Whether to apply VTST input swaps. Defaults to False in settings.
    vasp_custodian_handlers
        List of handlers to use in Custodian. See settings for list.
    vasp_custodian_validators
        List of validators to use in Custodian. See settings for list.
    scratch_dir
        Scratch directory to use. Defaults to None.
    vasp_job_kwargs
        Keyword arguments to pass to VaspJob. Defaults to None.
    custodian_kwargs
        Any remaining keyword arguments to pass to Custodian. Defaults to None.
    """
    # Adapted from https://github.com/materialsproject/atomate2/blob/main/src/atomate2/vasp/run.py

    vasp_parallel_cmd = os.path.expandvars(vasp_parallel_cmd)

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
        "ScanMetalHandler": ScanMetalHandler(),
    }
    validators_dict = {
        "VaspFilesValidator": VaspFilesValidator(),
        "VasprunXMLValidator": VasprunXMLValidator(),
    }

    handlers = []
    for handler_flag in vasp_custodian_handlers:
        if handler_flag not in handlers_dict:
            raise ValueError(f"Unknown VASP error handler: {handler_flag}")
        handlers.append(handlers_dict[handler_flag])

    validators = []
    for validator_flag in vasp_custodian_validators:
        if validator_flag not in validators_dict:
            raise ValueError(f"Unknown VASP validator: {validator_flag}")
        validators.append(validators_dict[validator_flag])

    # Populate settings
    vasp_cmd = f"{vasp_parallel_cmd} {vasp_cmd}"
    vasp_gamma_cmd = f"{vasp_parallel_cmd} {vasp_gamma_cmd}"

    # Run VASP
    vasp_job_kwargs = {} if vasp_job_kwargs is None else vasp_job_kwargs
    custodian_kwargs = {} if custodian_kwargs is None else custodian_kwargs
    split_vasp_cmd = shlex.split(vasp_cmd)
    split_vasp_gamma_cmd = shlex.split(vasp_gamma_cmd)
    vasp_job_kwargs.update({"gamma_vasp_cmd": split_vasp_gamma_cmd})

    # Run with Custodian
    jobs = [VaspJob(split_vasp_cmd, **vasp_job_kwargs)]

    if vasp_custodian_wall_time is not None:
        handlers = list(handlers) + [
            WalltimeHandler(wall_time=vasp_custodian_wall_time)
        ]

    c = Custodian(
        handlers,
        jobs,
        validators=validators,
        max_errors=vasp_custodian_max_errors,
        scratch_dir=scratch_dir,
        **custodian_kwargs,
    )

    c.run()


if __name__ == "__main__":
    run_custodian()
