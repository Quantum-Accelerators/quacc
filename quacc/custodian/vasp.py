import shlex

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


def run_custodian():
    # Adapted from https://github.com/materialsproject/atomate2/blob/main/src/atomate2/vasp/run.py

    # Handlers for VASP
    handlers = []
    handlers_dict = {
        "VaspErrorHandler": VaspErrorHandler(vtst_fixes=SETTINGS.VASP_CUSTODIAN_VTST),
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
    for handler_flag in SETTINGS.VASP_CUSTODIAN_HANDLERS:
        if handler_flag not in handlers_dict:
            raise ValueError(f"Unknown VASP error handler: {handler_flag}")
        handlers.append(handlers_dict[handler_flag])

    validators = []
    for validator_flag in SETTINGS.VASP_CUSTODIAN_VALIDATORS:
        if validator_flag not in validators_dict:
            raise ValueError(f"Unknown VASP validator: {validator_flag}")
        validators.append(validators_dict[validator_flag])

    # Populate settings
    parallel_cmd = SETTINGS.VASP_PARALLEL_CMD + " "
    vasp_cmd = parallel_cmd + SETTINGS.VASP_CMD
    vasp_gamma_cmd = parallel_cmd + SETTINGS.VASP_GAMMA_CMD
    max_errors = SETTINGS.VASP_CUSTODIAN_MAX_ERRORS
    wall_time = SETTINGS.VASP_CUSTODIAN_WALL_TIME
    scratch_dir = None
    vasp_job_kwargs = None
    custodian_kwargs = None

    # Run VASP
    vasp_job_kwargs = {} if vasp_job_kwargs is None else vasp_job_kwargs
    custodian_kwargs = {} if custodian_kwargs is None else custodian_kwargs
    split_vasp_cmd = shlex.split(vasp_cmd)
    split_vasp_gamma_cmd = shlex.split(vasp_gamma_cmd)
    vasp_job_kwargs.update({"gamma_vasp_cmd": split_vasp_gamma_cmd})

    # Run with Custodian
    jobs = [VaspJob(split_vasp_cmd, **vasp_job_kwargs)]

    if wall_time is not None:
        handlers = list(handlers) + [WalltimeHandler(wall_time=wall_time)]

    c = Custodian(
        handlers,
        jobs,
        validators=validators,
        max_errors=max_errors,
        scratch_dir=scratch_dir,
        **custodian_kwargs,
    )

    c.run()


if __name__ == "__main__":
    run_custodian()
