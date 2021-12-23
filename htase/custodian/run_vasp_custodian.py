import os
import shlex
import subprocess
import yaml
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

# Adapted from https://github.com/materialsproject/atomate2/blob/main/src/atomate2/vasp/run.py

# Read in default settings
if "VASP_CUSTODIAN_SETTINGS" in os.environ:
    settings_path = os.environ["VASP_CUSTODIAN_SETTINGS"]
else:
    settings_path = os.path.join(
        os.path.dirname(os.environ["RUN_VASP_CUSTODIAN"]),
        "vasp_custodian_settings.yaml",
    )
if not os.path.exists(settings_path):
    raise ValueError(
        "Missing vasp_custodian_settings.yaml in same directory as run_custodian.py"
    )
config = yaml.safe_load(open(settings_path))

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
custodian_enabled = config["custodian_enabled"]
vasp_cmd = config["vasp_cmd"]
vasp_gamma_cmd = config["vasp_gamma_cmd"]
max_errors = config["max_errors"]
scratch_dir = config["scratch_dir"]
wall_time = config["wall_time"]
vasp_job_kwargs = config["vasp_job_kwargs"]
custodian_kwargs = config["custodian_kwargs"]

vasp_job_kwargs = {} if vasp_job_kwargs is None else vasp_job_kwargs
custodian_kwargs = {} if custodian_kwargs is None else custodian_kwargs

vasp_cmd = os.path.expandvars(vasp_cmd)
vasp_gamma_cmd = os.path.expandvars(vasp_gamma_cmd)
split_vasp_cmd = shlex.split(vasp_cmd)
split_vasp_gamma_cmd = shlex.split(vasp_gamma_cmd)

vasp_job_kwargs["auto_npar"] = bool(vasp_job_kwargs.get("auto_npar", None))

vasp_job_kwargs.update({"gamma_vasp_cmd": split_vasp_gamma_cmd})

if custodian_enabled:

    # Run VASP with custodian
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

else:

    # Run VASP without custodian
    return_code = subprocess.call(vasp_cmd, shell=True)
