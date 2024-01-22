"""Settings for quacc."""
from __future__ import annotations

import os
from importlib import util
from pathlib import Path
from shutil import which
from typing import TYPE_CHECKING, Literal, Optional, Union

import psutil
from maggma.core import Store
from pydantic import Field, field_validator, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

if TYPE_CHECKING:
    from typing import Any

installed_engine = next(
    (
        wflow_engine
        for wflow_engine in ["parsl", "covalent", "dask", "redun", "jobflow"]
        if util.find_spec(wflow_engine)
    ),
    None,
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
    using the "QUACC" prefix. e.g. `export QUACC_SCRATCH_DIR=/path/to/scratch`.
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

    WORKFLOW_ENGINE: Optional[
        Literal["covalent", "dask", "parsl", "prefect", "redun", "jobflow"]
    ] = Field(installed_engine, description=("The workflow manager to use, if any."))

    # ---------------------------
    # General Settings
    # ---------------------------

    RESULTS_DIR: Path = Field(
        Path.cwd(),
        description=(
            "Directory to permanently store I/O-based calculation results in."
            "Note that this behavior may be modified by the chosen workflow engine."
        ),
    )
    SCRATCH_DIR: Optional[Path] = Field(
        None,
        description=(
            "The base directory where calculations are run. If set to None, calculations will be run in a "
            "temporary directory within `RESULTS_DIR`. If a `Path` is supplied, calculations will "
            "be run in a temporary directory within `SCRATCH_DIR`. Files are always moved back "
            "to `RESULTS_DIR` after the calculation is complete, and the temporary directory "
            "in `SCRATCH_DIR` is removed."
        ),
    )
    CREATE_UNIQUE_DIR: bool = Field(
        True,
        description=(
            "Whether to have a unique directory in RESULTS_DIR for each job."
            "Some workflow engines have an option to do this for you already,"
            "in which case you should set this to False."
        ),
    )
    GZIP_FILES: bool = Field(
        True, description="Whether generated files should be gzip'd."
    )
    CHECK_CONVERGENCE: bool = Field(
        True,
        description="Whether to check for convergence, when implemented by a given recipe.",
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
    # Prefect Settings
    # ---------------------------
    PREFECT_AUTO_SUBMIT: bool = Field(
        True, description="Whether to auto-submit tasks to the task runner."
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
    # ESPRESSO Settings
    # ---------------------------
    ESPRESSO_BIN_DIR: Path = Field(
        Path(), description="Base path to the espresso binaries."
    )
    ESPRESSO_BINARIES: dict[str, str] = Field(
        {
            "pw": "pw.x",
            "ph": "ph.x",
            "neb": "neb.x",
            "q2r": "q2r.x",
            "dos": "dos.x",
            "matdyn": "matdyn.x",
            "dynmat": "dynmat.x",
            "bands": "bands.x",
            "projwfc": "projwfc.x",
            "pp": "pp.x",
            "wannier90": "wannier90.x",
        },
        description="Name for each espresso binary.",
    )
    ESPRESSO_PSEUDO: Optional[Path] = Field(
        None, description=("Path to a pseudopotential library for espresso.")
    )
    ESPRESSO_PRESET_DIR: Path = Field(
        Path(__file__).parent / "calculators" / "espresso" / "presets",
        description="Path to the espresso preset directory",
    )

    # ---------------------------
    # Gaussian Settings
    # ---------------------------
    GAUSSIAN_CMD: Path = Field(
        Path("g16"), description=("Path to the Gaussian executable.")
    )

    # ---------------------------
    # ONETEP Settings
    # ---------------------------
    ONETEP_CMD: Optional[Path] = Field(
        Path("onetep.arch"), description=("Path to the ONETEP executable.")
    )
    ONETEP_PARALLEL_CMD: Optional[dict] = Field(
        None,
        description=(
            "Parallelization commands to run ONETEP that are prepended to the executable."
        ),
    )
    ONETEP_PP_PATH: Optional[Path] = Field(
        None, description=("Path to pseudopotentials.")
    )

    # ---------------------------
    # GULP Settings
    # ---------------------------
    GULP_CMD: Path = Field(Path("gulp"), description=("Path to the GULP executable."))
    GULP_LIB: Optional[Path] = Field(
        None,
        description=(
            "Path to the GULP force field library. If not specified, the GULP_LIB environment variable will be used (if present)."
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
    QCHEM_NUM_CORES: int = Field(
        psutil.cpu_count(logical=False),
        description="Number of cores to use for the Q-Chem calculation.",
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

    # ---------------------------
    # Debug Settings
    # ---------------------------
    DEBUG: bool = Field(
        False,
        description=(
            "Whether to run in debug mode. This will set the logging level to DEBUG, "
            "ASE logs (e.g. optimizations, vibrations, thermo) are printed to stdout."
        ),
    )

    # --8<-- [end:settings]

    @field_validator("RESULTS_DIR", "SCRATCH_DIR")
    @classmethod
    def resolve_and_make_paths(cls, v: Optional[Path]) -> Optional[Path]:
        """Resolve and make paths."""
        if v is None:
            return v

        v = Path(os.path.expandvars(v)).expanduser().resolve()
        if not v.exists():
            v.mkdir(parents=True)
        return v

    @field_validator(
        "ESPRESSO_PRESET_DIR",
        "ESPRESSO_PSEUDO",
        "GAUSSIAN_CMD",
        "GULP_CMD",
        "GULP_LIB",
        "ORCA_CMD",
        "QCHEM_LOCAL_SCRATCH",
        "NEWTONNET_MODEL_PATH",
        "VASP_PRESET_DIR",
        "VASP_PP_PATH",
        "VASP_VDW",
    )
    @classmethod
    def expand_paths(cls, v: Optional[Path]) -> Optional[Path]:
        """Expand ~/ in paths."""
        return v.expanduser() if v is not None else v

    @field_validator("PRIMARY_STORE")
    def generate_store(cls, v: Union[str, Store]) -> Store:
        """Generate the Maggma store."""
        from monty.json import MontyDecoder

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
