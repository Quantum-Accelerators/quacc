"""Settings for quacc"""
from __future__ import annotations

import os
from shutil import which
from typing import List, Literal, Optional

from pydantic import BaseSettings, Field, root_validator

from quacc.presets import vasp as vasp_defaults

_DEFAULT_CONFIG_FILE_PATH = os.path.join(os.path.expanduser("~"), ".quacc.yaml")


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

    WORKFLOW_MANAGER: Optional[Literal["covalent", "parsl"]] = Field(
        "covalent", description="The workflow manager to use."
    )
    CONFIG_FILE: str = Field(
        _DEFAULT_CONFIG_FILE_PATH, description="File to load alternative defaults from."
    )
    RESULTS_DIR: str = Field(
        os.getcwd(),
        description="Directory to store results in.",
    )
    SCRATCH_DIR: str = Field(
        os.path.join("/tmp", "quacc") if os.path.exists("/tmp") else os.getcwd(),
        description="Scratch directory for calculations.",
    )
    CREATE_UNIQUE_WORKDIR: bool = Field(
        False,
        description="Whether to automatically create a unique working directory for each calculation. Some workflow engines have an option to do this for you already.",
    )
    GZIP_FILES: bool = Field(
        True, description="Whether generated files should be gzip'd."
    )

    # ---------------------------
    # Data Store Settings
    # ---------------------------
    PRIMARY_STORE: str = Field(
        None,
        description="String-based JSON representation of the primary Maggma data store where calculation results will be stored. Taken from the `.to_json()` method of the corresponding Store object.",
    )

    # ---------------------------
    # ORCA Settings
    # ---------------------------
    ORCA_CMD: str = Field(
        "orca",
        description="Path to the ORCA executable. This must be the full, absolute path for parallel calculations to work.",
    )

    # ---------------------------
    # VASP Settings
    # ---------------------------

    # VASP Settings: Main
    VASP_PARALLEL_CMD: str = Field(
        "",
        description="Parallel command to run VASP with Custodian (e.g. srun -N 2 --ntasks-per-node 48)",
    )
    VASP_CMD: str = Field(
        "vasp_std", description="Command to run the standard version of VASP."
    )
    VASP_GAMMA_CMD: str = Field(
        "vasp_gam", description="Command to run the gamma-point only version of VASP."
    )

    # VASP Settings: General
    VASP_INCAR_COPILOT: bool = Field(
        True, description="Whether co-pilot mode should be used for VASP INCAR handling"
    )
    VASP_BADER: bool = Field(
        bool(which("bader")),
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
    VASP_COPY_MAGMOMS: bool = Field(
        True,
        description="If True, any pre-existing atoms.get_magnetic_moments() will be set in atoms.set_initial_magnetic_moments().",
    )
    VASP_VERBOSE: bool = Field(
        True,
        description="If True, warnings will be raised when INCAR parameters are changed.",
    )
    VASP_PRESET_DIR: str = Field(
        os.path.dirname(vasp_defaults.__file__),
        description="Path to the VASP preset directory",
    )

    # VASP Settings: Custodian
    VASP_USE_CUSTODIAN: bool = Field(
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
    VASP_CUSTODIAN_WALL_TIME: Optional[int] = Field(
        None,
        description="After this many seconds, Custodian will stop running and ensure that VASP writes a STOPCAR",
    )

    class Config:
        """Pydantic config settings."""

        env_prefix = "quacc_"

    @root_validator(pre=True)
    def load_default_settings(cls, values: dict) -> dict:
        """
        Loads settings from a root file if available and uses that as defaults in
        place of built in defaults.

        Parameters
        ----------
        values
            Settings to load.

        Returns
        -------
        dict
            Loaded settings.
        """

        from monty.serialization import loadfn

        config_file_path = values.get("CONFIG_FILE", _DEFAULT_CONFIG_FILE_PATH)

        new_values = {}
        if os.path.exists(os.path.expanduser(config_file_path)):
            new_values |= loadfn(os.path.expanduser(config_file_path))

        new_values.update(values)
        return new_values
