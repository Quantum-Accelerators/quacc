import os
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
from htase.util.custodian import load_yaml_settings

# Adapted from https://github.com/materialsproject/atomate2/blob/main/src/atomate2/vasp/run.py

# Read in default settings
if "VASP_CUSTODIAN_SETTINGS" in os.environ:
    settings_path = os.environ["VASP_CUSTODIAN_SETTINGS"]
else:
    raise EnvironmentError("Missing environment variable VASP_CUSTODIAN_SETTINGS.")
config = load_yaml_settings(settings_path)

# Handlers for VASP
handlers = []
handlers_dict = {
    "VaspErrorHandler": VaspErrorHandler(vtst_fixes=config["vtst_swaps"]),
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
for handler_flag in config["handlers"]:
    if handler_flag not in handlers_dict.keys():
        raise ValueError(f"Unknown VASP error handler: {handler_flag}")
    handlers.append(handlers_dict[handler_flag])

validators = []
for validator_flag in config["validators"]:
    if validator_flag not in validators_dict.keys():
        raise ValueError(f"Unknown VASP validator: {validator_flag}")
    validators.append(validators_dict[validator_flag])

# Populate settings
parallel_cmd = config.get("vasp_parallel_cmd", "") + " "
vasp_cmd = parallel_cmd + config.get("vasp_cmd", "vasp_std")
vasp_gamma_cmd = parallel_cmd + config.get("vasp_gamma_cmd", "vasp_gam")
max_errors = config.get("max_errors", 5)
wall_time = config.get("custodian_wall_time", None)
scratch_dir = config.get("scratch_dir", None)
vasp_job_kwargs = config.get("vasp_job_kwargs", None)
custodian_kwargs = config.get("custodian_kwargs", None)

# Run VASP
vasp_job_kwargs = {} if vasp_job_kwargs is None else vasp_job_kwargs
custodian_kwargs = {} if custodian_kwargs is None else custodian_kwargs
vasp_cmd = os.path.expandvars(vasp_cmd)
vasp_gamma_cmd = os.path.expandvars(vasp_gamma_cmd)
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

print("Running VASP using custodian.")
c.run()
