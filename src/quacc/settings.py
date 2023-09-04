"""Settings for quacc"""
from __future__ import annotations

import os
from pathlib import Path
from shutil import which
from typing import List, Optional, Union

from pydantic import BaseSettings, Field, root_validator

from quacc.presets import vasp as vasp_defaults

try:
    import covalent
except ImportError:
    covalent = None
try:
    import parsl
except ImportError:
    parsl = None
try:
    import jobflow
except ImportError:
    jobflow = None
try:
    import redun
except ImportError:
    redun = None

_DEFAULT_CONFIG_FILE_PATH = Path("~", ".quacc.yaml").expanduser()


class QuaccSettings(BaseSettings):
    """
    Settings for quacc.

    The default way to modify these is to make a ~/.quacc.yaml file.
    Alternatively, the environment variable QUACC_CONFIG_FILE can be set to
    point to a yaml file with quacc settings.

    The variables can also be modified individually though environment variables
    by using the "QUACC" prefix. e.g. QUACC_SCRATCH_DIR=/path/to/scratch.
    """

    # --8<-- [start:settings]

    # ---------------------------
    # Workflow Engine
    # ---------------------------

    WORKFLOW_ENGINE: str = Field(
        "covalent"
        if covalent
        else "parsl"
        if parsl
        else "redun"
        if redun
        else "jobflow"
        if jobflow
        else "local",
        description=(
            "The workflow manager to use."
            "Options include: 'covalent', 'parsl', 'redun', 'jobflow', or 'local'"
        ),
    )

    # ---------------------------
    # General Settings
    # ---------------------------

    CONFIG_FILE: Union[str, Path] = Field(
        _DEFAULT_CONFIG_FILE_PATH,
        description=(
            "Path to the YAML file to load alternative quacc configuration "
            "defaults from."
        ),
    )
    RESULTS_DIR: Union[str, Path] = Field(
        Path.cwd(),
        description=(
            "Directory to store I/O-based calculation results in."
            "Note that this behavior may be modified by the chosen workflow engine."
            "For instance, Covalent specifies the base directory as the `workdir` "
            "of a local executor or the `remote_workdir` of a remote executor."
            "In this case, the `RESULTS_DIR` will be a subdirectory of that directory."
        ),
    )
    SCRATCH_DIR: Union[str, Path] = Field(
        Path("/tmp") if Path("/tmp").exists() else Path.cwd(),
        description="Scratch directory for calculations.",
    )
    CREATE_UNIQUE_WORKDIR: bool = Field(
        False,
        description=(
            "Whether to have a unique working directory in RESULTS_DIR for each job."
            "Some workflow engines have an option to do this for you already."
        ),
    )
    GZIP_FILES: bool = Field(
        True, description="Whether generated files should be gzip'd."
    )
    CHECK_CONVERGENCE: bool = Field(
        True,
        description="Whether to check for convergence in the `summarize_run`-type functions, if supported.",
    )

    # ---------------------------
    # Data Store Settings
    # ---------------------------
    PRIMARY_STORE: str = Field(
        None,
        description=(
            "String-based JSON representation of the primary Maggma data store "
            "where calculation results will be stored."
            "Taken from the `.to_json()` method of the corresponding Store object."
        ),
    )

    # ---------------------------
    # ORCA Settings
    # ---------------------------
    ORCA_CMD: Union[str, Path] = Field(
        "orca",
        description=(
            "Path to the ORCA executable. This must be the full, absolute path "
            "for parallel calculations to work."
        ),
    )

    # ---------------------------
    # VASP Settings
    # ---------------------------

    # VASP Settings: Main
    VASP_PARALLEL_CMD: str = Field(
        "",
        description=(
            "Parallel command to run VASP with Custodian."
            "For example: srun -N 2 --ntasks-per-node 48"
            "Note that this does not include the executable name."
        ),
    )
    VASP_CMD: str = Field(
        "vasp_std", description="Command to run the standard version of VASP."
    )
    VASP_GAMMA_CMD: str = Field(
        "vasp_gam", description="Command to run the gamma-point only version of VASP."
    )

    # VASP Settings: General
    VASP_INCAR_COPILOT: bool = Field(
        True,
        description=(
            "Whether co-pilot mode should be used for VASP INCAR handling."
            "This will modify INCAR flags on-the-fly if they disobey the VASP manual."
            "A warning will be raised in each case."
        ),
    )
    VASP_BADER: bool = Field(
        bool(which("bader")),
        description=(
            "Whether to run a Bader analysis when summarizing VASP results."
            "Requires bader to be in PATH."
        ),
    )
    VASP_PRESET_MAG_DEFAULT: float = Field(
        1.0,
        description=(
            "Default initial magmom to use for a given element if a preset "
            "with magmoms is provided but an element is missing from the list"
        ),
    )
    VASP_MAG_CUTOFF: float = Field(
        0.05,
        description=(
            "If the absolute value of all magnetic moments are below this value, "
            "they will be set to 0 such that a spin-unpolarized calculation will be performed"
        ),
    )
    VASP_COPY_MAGMOMS: bool = Field(
        True,
        description=(
            "If True, any pre-existing atoms.get_magnetic_moments() will be set"
            "in atoms.set_initial_magnetic_moments()."
        ),
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
        description=(
            "If VTST-related input swaps should be used when running Custodian."
            "Requires VASP to be compiled with VTST"
        ),
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
        description=(
            "After this many seconds, Custodian will stop running "
            "and ensure that VASP writes a STOPCAR"
        ),
    )

    # ---------------------------
    # Q-Chem Settings
    # ---------------------------

    # Q-Chem Settings: Main
    QCHEM_CMD: str = Field(
        "qchem", description="Command to run the standard version of Q-Chem."
    )

    QCHEM_LOCAL_SCRATCH: Union[str, Path] = Field(
        Path("/tmp") if Path("/tmp").exists() else Path.cwd(),
        description="Compute-node local scratch directory in which Q-Chem should perform IO.",
    )

    # Q-Chem Settings: Custodian
    QCHEM_USE_ERROR_HANDLERS: bool = Field(
        True,
        description="Whether Custodian's error handlers should be employed for Q-Chem.",
    )

    QCHEM_CUSTODIAN_MAX_ERRORS: int = Field(
        5, description="Maximum errors for Q-Chem Custodian."
    )

    # ---------------------------
    # NewtonNet Settings
    # ---------------------------
    NEWTONNET_MODEL_PATH: Union[Union[str, Path], List[Union[str, Path]]] = Field(
        "best_model_state.tar", description="Path to NewtonNet .tar model"
    )
    NEWTONNET_CONFIG_PATH: Union[Union[str, Path], List[Union[str, Path]]] = Field(
        "config.yml", description="Path to NewtonNet YAML settings file"
    )

    # --8<-- [end:settings]

    class Config:
        """Pydantic config settings."""

        env_prefix = "quacc_"

    @root_validator(pre=True)
    def load_default_settings(cls, values: dict) -> dict:
        """
        Loads settings from a root file if available and uses that as defaults
        in place of built in defaults.

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

        config_file_path = Path(
            values.get("CONFIG_FILE", _DEFAULT_CONFIG_FILE_PATH)
        ).expanduser()

        new_values = {}
        if config_file_path.exists() and config_file_path.stat().st_size > 0:
            new_values |= loadfn(config_file_path)

        new_values.update(values)
        return new_values
