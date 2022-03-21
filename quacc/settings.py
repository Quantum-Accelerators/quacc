"""Settings for quacc"""
import os
from pathlib import Path
from typing import List
from pydantic import BaseSettings, Field, root_validator


_DEFAULT_CONFIG_FILE_PATH = "~/.quacc.yaml"


class QuaccSettings(BaseSettings):
    """
    Settings for quacc.

    The default way to modify these is to make a ~/.quacc.yaml file. Alternatively,
    the environment variable QUACC_CONFIG_FILE can be set to point to a yaml file with
    quacc settings.

    The variables can also be modified individually though environment variables by
    using the "QUACC" prefix. e.g. QUACC_SCRATCH_DIR=/path/to/scratch.
    """

    # ---------------------------
    # General Settings
    # ---------------------------

    CONFIG_FILE: str = Field(
        _DEFAULT_CONFIG_FILE_PATH, description="File to load alternative defaults from."
    )
    SCRATCH_DIR: str = Field(
        os.path.expandvars("$SCRATCH") if "SCRATCH" in os.environ else ".",
        description="Scratch directory for calculations.",
    )
    GZIP_FILES: bool = Field(
        True, description="Whether generated files should be gzip'd."
    )

    # ---------------------------
    # VASP Settings
    # ---------------------------

    # VASP Settings: Main
    VASP_PARALLEL_CMD: str = Field(
        os.path.expandvars("$VASP_PARALLEL_CMD"),
        description="Parallel command to run VASP with Custodian (e.g. srun -N 2 --ntasks-per-node 24)",
    )
    VASP_CMD: str = Field(
        "vasp_std", description="Command to run the standard version of VASP."
    )
    VASP_GAMMA_CMD: str = Field(
        "vasp_gam", description="Command to run the gamma-point only version of VASP."
    )

    # VASP Settings: General
    INCAR_COPILOT: bool = Field(
        True, description="Whether co-pilot mode should be used for VASP INCAR handling"
    )
    VASP_BADER: bool = Field(
        True,
        description="Whether to run a Bader analysis when summarizing VASP results. Requires bader to be in PATH.",
    )
    VASP_PRESET_MAG_DEFAULT: float = Field(
        1.0,
        description="Default initial magmom to use for a given element if a preset with magmoms is provided but an element is missing from the list",
    )
    VASP_MAG_CUTOFF: float = Field(
        0.05,
        description="If the absolute value of all magnetic moments are below this value, they will be set to 0 such that a spin-unpolarized calculation will be performed",
    )

    # VASP Settings: Custodian
    VASP_CUSTODIAN: bool = Field(
        True, description="Whether Custodian should be used to run VASP"
    )
    VASP_CUSTODIAN_VTST: bool = Field(
        False,
        description="If VTST-related input swaps should be used when running Custodian. Requires VASP to be compiled with VTST",
    )
    VASP_CUSTODIAN_MAX_ERRORS: int = Field(
        5, description="Maximum errors for Custodian"
    )
    VASP_CUSTODIAN_HANDLERS: List[str] = Field(
        [
            "VaspErrorHandler",
            "MeshSymmetryErrorHandler",
            "UnconvergedErrorHandler",
            "NonConvergingErrorHandler",
            "PotimErrorHandler",
            "PositiveEnergyErrorHandler",
            "FrozenJobErrorHandler",
            "StdErrHandler",
            "LargeSigmaHandler",
            "IncorrectSmearingHandler",
        ],
        description="Handlers for Custodian",
    )
    VASP_CUSTODIAN_VALIDATORS: List[str] = Field(
        ["VasprunXMLValidator", "VaspFilesValidator"],
        description="Validators for Custodian",
    )
    VASP_CUSTODIAN_WALL_TIME: int = Field(
        None,
        description="After this many seconds, Custodian will stop running and ensure that VASP writes a STOPCAR",
    )

    class Config:
        """Pydantic config settings."""

        env_prefix = "quacc_"

    @root_validator(pre=True)
    def load_default_settings(cls, values):
        """
        Load settings from file or environment variables.
        Loads settings from a root file if available and uses that as defaults in
        place of built in defaults.
        This allows setting of the config file path through environment variables.
        """
        from monty.serialization import loadfn

        config_file_path: str = values.get("CONFIG_FILE", _DEFAULT_CONFIG_FILE_PATH)

        new_values = {}
        if Path(config_file_path).exists():
            new_values.update(loadfn(config_file_path))

        new_values.update(values)
        return new_values
