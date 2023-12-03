"""Settings for quacc."""
from __future__ import annotations

import os
from importlib import util
from pathlib import Path
from shutil import which
from typing import TYPE_CHECKING, Literal, Optional, Union

from maggma.core import Store
from monty.json import MontyDecoder
from pydantic import Field, field_validator, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

if TYPE_CHECKING:
    from typing import Any

installed_engine = next(
    (
        wflow_engine
        for wflow_engine in ["parsl", "covalent", "prefect", "redun", "jobflow"]
        if util.find_spec(wflow_engine)
    ),
    "local",
)
_DEFAULT_CONFIG_FILE_PATH = Path("~", ".quacc.yaml").expanduser().resolve()


class QuaccSettings(BaseSettings):
    """
    Settings for quacc.

    The default way to modify these is to make a ~/.quacc.yaml file. Alternatively, the
    environment variable QUACC_CONFIG_FILE can be set to point to a custom yaml file
    with quacc settings. The quacc CLI offers a `quacc set <setting> <value>` option to
    do this as well.

    The variables can also be modified individually though environment variables by
    using the "QUACC" prefix. e.g. QUACC_SCRATCH_DIR=/path/to/scratch.
    """

    CONFIG_FILE: Path = Field(
        _DEFAULT_CONFIG_FILE_PATH,
        description=(
            "Path to the YAML file to load alternative quacc configuration "
            "defaults from."
        ),
    )

    # --8<-- [start:settings]

    # ---------------------------
    # Workflow Engine
    # ---------------------------

    WORKFLOW_ENGINE: Literal[
        "covalent", "jobflow", "parsl", "prefect", "redun", "local"
    ] = Field(
        installed_engine,
        description=(
            "The workflow manager to use."
            "Options include: 'covalent', 'parsl', 'redun', 'jobflow', 'prefect', or 'local'"
        ),
    )

    # ---------------------------
    # General Settings
    # ---------------------------

    RESULTS_DIR: Path = Field(
        Path.cwd(),
        description=(
            "Directory to store I/O-based calculation results in."
            "Note that this behavior may be modified by the chosen workflow engine."
            "For instance, Covalent specifies the base directory as the `workdir` "
            "of a local executor or the `remote_workdir` of a remote executor."
            "In this case, the `RESULTS_DIR` will be a subdirectory of that directory."
        ),
    )
    SCRATCH_DIR: Path = Field(
        Path("~/.scratch"), description="Scratch directory for calculations."
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
    PRIMARY_STORE: Optional[Union[str, Store]] = Field(
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
    ORCA_CMD: Path = Field(
        Path(which("orca") or "orca"),
        description=(
            "Path to the ORCA executable. This must be the full, absolute path "
            "for parallel calculations to work."
        ),
    )

    # ---------------------------
    # Gaussian Settings
    # ---------------------------
    GAUSSIAN_CMD: Path = Field(
        Path("g16"), description=("Path to the Gaussian executable.")
    )

    # ---------------------------
    # GULP Settings
    # ---------------------------
    GULP_CMD: Path = Field(Path("gulp"), description=("Path to the GULP executable."))
    GULP_LIB: Optional[Path] = Field(
        None, description=("Path to the GULP force field library.")
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
    VASP_PP_PATH: Optional[Path] = Field(
        None,
        description="Path to the VASP pseudopotential library. Must contain the directories `potpaw_PBE` and `potpaw` for PBE and LDA pseudopotentials, respectively.",
    )
    VASP_VDW: Optional[Path] = Field(
        None, description="Path to the vdw_kernel.bindat file for VASP vdW functionals."
    )

    # VASP Settings: General
    VASP_INCAR_COPILOT: Literal["off", "on", "aggressive"] = Field(
        "on",
        description=(
            "Controls VASP co-pilot mode for automated INCAR parameter handling."
            "off: Do not use co-pilot mode. INCAR parameters will be unmodified."
            "on: Use co-pilot mode. This will only modify INCAR flags not already set by the user."
            "aggressive: Use co-pilot mode in aggressive mode. This will modify INCAR flags even if they are already set by the user."
        ),
    )
    VASP_BADER: bool = Field(
        bool(which("bader")),
        description=(
            "Whether to run a Bader analysis when summarizing VASP results."
            "Requires bader to be in PATH."
        ),
    )
    VASP_CHARGEMOL: bool = Field(
        bool(os.environ.get("DDEC6_ATOMIC_DENSITIES_DIR")),
        description=(
            "Whether to run a Chargemol (i.e. DDEC6, CM5) analysis when summarizing VASP results."
            "Requires the Chargemol executable to be in PATH and the DDEC6_ATOMIC_DENSITIES_DIR environment variable."
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
    VASP_PRESET_DIR: Path = Field(
        Path(__file__).parent / "calculators" / "vasp" / "presets",
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
    VASP_CUSTODIAN_HANDLERS: list[str] = Field(
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
    VASP_CUSTODIAN_VALIDATORS: list[str] = Field(
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

    QCHEM_LOCAL_SCRATCH: Path = Field(
        Path("/tmp") if Path("/tmp").exists() else Path.cwd() / ".qchem_scratch",
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

    # NBO Settings
    QCHEM_NBO_EXE: Optional[Path] = Field(
        None, description="Full path to the NBO executable."
    )

    # ---------------------------
    # NewtonNet Settings
    # ---------------------------
    NEWTONNET_MODEL_PATH: Union[Path, list[Path]] = Field(
        "best_model_state.tar", description="Path to NewtonNet .tar model"
    )
    NEWTONNET_CONFIG_PATH: Union[Path, list[Path]] = Field(
        "config.yml", description="Path to NewtonNet YAML settings file"
    )

    # --8<-- [end:settings]

    @field_validator("RESULTS_DIR", "SCRATCH_DIR")
    @classmethod
    def resolve_and_make_paths(cls, v):
        v = Path(os.path.expandvars(v)).expanduser().resolve()
        if not v.exists():
            os.makedirs(v)
        return v

    @field_validator(
        "GAUSSIAN_CMD", "ORCA_CMD", "QCHEM_LOCAL_SCRATCH", "VASP_PRESET_DIR"
    )
    @classmethod
    def expand_paths(cls, v):
        return v.expanduser()

    @field_validator("PRIMARY_STORE")
    def generate_store(cls, v):
        return MontyDecoder().decode(v) if isinstance(v, str) else v

    model_config = SettingsConfigDict(env_prefix="quacc_")

    @model_validator(mode="before")
    @classmethod
    def load_default_settings(cls, values: dict[str, Any]) -> dict[str, Any]:
        """
        Loads settings from a root file if available and uses that as defaults in place
        of built in defaults.

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

        config_file_path = (
            Path(values.get("CONFIG_FILE", _DEFAULT_CONFIG_FILE_PATH))
            .expanduser()
            .resolve()
        )

        new_values = {}  # type: dict
        if config_file_path.exists() and config_file_path.stat().st_size > 0:
            new_values |= loadfn(config_file_path)

        new_values.update(values)
        return new_values
