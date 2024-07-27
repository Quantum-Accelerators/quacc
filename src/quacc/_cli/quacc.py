"""Quacc CLI module."""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from rich import print as rich_print
from typer import Exit, Option, Typer

from quacc import get_settings
from quacc.settings import QuaccSettings, _type_handler

app = Typer()

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Any


def callback(value: bool) -> None:
    """
    Set up the callback for the quacc version reporting.

    Parameters
    ----------
    value
        If the version should be reported

    Returns
    -------
    None
    """
    from quacc import __version__

    if value:
        rich_print(f"quacc v{__version__}")
        raise Exit


@app.callback()
def main(
    version: Optional[bool] = Option(  # noqa: ARG001, UP007
        None,
        "--version",
        "-v",
        help="Show the application's version and exit.",
        callback=callback,
        is_eager=True,
    ),
) -> None:
    """
    The main CLI interface, with an option to return the version.

    Parameters
    ----------
    version
        If the version should be reported

    Returns
    -------
    None
    """


@app.command("set")
def set_(parameter: str, new_value: str) -> None:
    """
    Set the specified quacc parameter in the quacc configuration file. This command will
    not override any environment variables.

    Parameters
    ----------
    parameter
        The quacc parameter to set.
    new_value
        The value of the quacc parameter.

    Returns
    -------
    None
    """
    from quacc.settings import _DEFAULT_CONFIG_FILE_PATH

    settings = get_settings()

    config_file = settings.CONFIG_FILE or _DEFAULT_CONFIG_FILE_PATH
    parameter = parameter.upper()

    new_value = _type_handler({parameter: new_value})[parameter]
    _parameter_handler(parameter, settings.model_dump(), value=new_value)

    rich_print(f"Setting `{parameter}` to `{new_value}` in {config_file}")
    _update_setting(parameter, new_value, config_file)


@app.command()
def unset(parameter: str) -> None:
    """
    Unset the specified quacc parameter in the quacc configuration file. This command
    will not override any environment variables.

    Parameters
    ----------
    parameter
        The quacc parameter to unset.

    Returns
    -------
    None
    """
    from quacc.settings import _DEFAULT_CONFIG_FILE_PATH

    settings = get_settings()
    CONFIG_FILE = settings.CONFIG_FILE or _DEFAULT_CONFIG_FILE_PATH
    parameter = parameter.upper()

    _parameter_handler(parameter, settings.model_dump())
    rich_print(f"Unsetting `{parameter}` in {CONFIG_FILE}")
    _delete_setting(parameter, CONFIG_FILE)


@app.command("info")
def info() -> None:
    """
    Print out some basic information about the quacc environment.

    Returns
    -------
    None
    """
    import platform

    from quacc import __version__

    settings = get_settings().model_dump()
    rich_print(
        f"""
quacc version: {__version__}

Python version: {platform.python_version()}

quacc settings: {settings}
"""
    )


def _parameter_handler(
    parameter: str, settings_dict: dict, value: Any | None = object
) -> None:
    """
    Check if the parameter is a valid quacc configuration variable.

    Parameters
    ----------
    parameter
        The quacc parameter to check.
    settings_dict
        The quacc settings.
    value
        The value to set the parameter to.

    Returns
    -------
    None
    """
    if parameter not in settings_dict:
        msg = f"{parameter} is not a supported quacc configuration variable."
        raise ValueError(msg)
    if parameter == "CONFIG_FILE":
        msg = "Cannot set the CONFIG_FILE parameter via the CLI."
        raise ValueError(msg)
    if value is not object:
        QuaccSettings(**{f"{parameter}": value})


def _delete_setting(key: str, config_file: Path) -> None:
    """
    Remove the quacc setting from the configuration file.

    Parameters
    ----------
    key
        The key in the YAML file to unset.
    config_file
        The path to the configuration file.

    Returns
    -------
    None
    """
    from ruamel.yaml import YAML

    yaml = YAML()
    if config_file.exists():
        with config_file.open() as yaml_file:
            yaml_content = yaml.load(yaml_file)

        if yaml_content:
            yaml_content.pop(key, None)
            with config_file.open(mode="w") as yaml_file:
                if yaml_content:
                    yaml.dump(yaml_content, yaml_file)


def _update_setting(key: str, value: Any, config_file: Path) -> None:
    """
    Update the quacc setting from the configuration file.

    Parameters
    ----------
    key
        The key in the YAML file to (re)set.
    value
        The value to set.
    config_file
        The path to the configuration file.

    Returns
    -------
    None
    """
    import ruamel.yaml

    yaml = ruamel.yaml.YAML()

    if config_file.exists() and config_file.stat().st_size > 0:
        with config_file.open() as yaml_file:
            yaml_content = yaml.load(yaml_file)
    else:
        yaml_content = {}

    yaml_content[key] = value
    with config_file.open(mode="w") as yaml_file:
        yaml.dump(yaml_content, yaml_file)


if __name__ == "__main__":
    app()
