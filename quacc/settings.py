"""Settings for quacc"""
import os
from pathlib import Path

from pydantic import BaseSettings, Field, root_validator

_DEFAULT_CONFIG_FILE_PATH = "~/.quacc.yaml"


class QuaccSettings(BaseSettings):
    """
    Settings for quacc.

    The default way to modify these is to modify ~/.quacc.yaml. Alternatively,
    the environment variable QUACC_CONFIG_FILE can be set to point to a yaml file with
    quacc settings.

    Alternatively, the variables can be modified directly though environment variables by
    using the "QUACC" prefix. e.g. QUACC_SCRATCH_DIR = /path/to/scratch.
    """

    CONFIG_FILE: str = Field(
        _DEFAULT_CONFIG_FILE_PATH, description="File to load alternative defaults from."
    )
    SCRATCH_DIR: str = Field(
        os.path.expandvars("$SCRATCH"),
        description="Scratch directory for calculations.",
    )
    GZIP_FILES: bool = Field(
        True, description="Whether generated files should be gzip'd."
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
