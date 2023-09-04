"""Quacc CLI module."""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from typing import Any

from pathlib import Path

import ruamel.yaml
import typer

from quacc import __version__
from quacc.settings import _DEFAULT_CONFIG_FILE_PATH

app = typer.Typer()


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
    if value:
        typer.echo(f"quacc v{__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        "-v",
        help="Show the application's version and exit.",
        callback=callback,
        is_eager=True,
    )
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
def set_(parameter: str, new_value) -> None:
    """
    Set the specified quacc parameter in the quacc configuration file. This
    command will not override any environment variables.

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
    from quacc import SETTINGS

    CONFIG_FILE = Path(SETTINGS.CONFIG_FILE or _DEFAULT_CONFIG_FILE_PATH)
    parameter = parameter.upper()
    if parameter not in SETTINGS.dict():
        msg = f"{parameter} is not a supported quacc configuration variable."
        raise ValueError(msg)
    if parameter == "CONFIG_FILE":
        msg = "Cannot set the CONFIG_FILE parameter via the CLI."
        raise ValueError(msg)

    typer.echo(f"Setting `{parameter}` to `{new_value}` in {CONFIG_FILE}")
    _update_setting(parameter, new_value, CONFIG_FILE)


@app.command()
def unset(parameter: str) -> None:
    """
    Unset the specified quacc parameter in the quacc configuration file. This
    command will not override any environment variables.

    Parameters
    ---------
    parameter
        The quacc parameter to unset.

    Returns
    -------
    None
    """
    from quacc import SETTINGS

    CONFIG_FILE = Path(SETTINGS.CONFIG_FILE or _DEFAULT_CONFIG_FILE_PATH)
    parameter = parameter.upper()
    if parameter not in SETTINGS.dict():
        msg = f"{parameter} is not a supported quacc configuration variable."
        raise ValueError(msg)
    if parameter == "CONFIG_FILE":
        msg = "Cannot unset the CONFIG_FILE parameter via the CLI."
        raise ValueError(msg)

    typer.echo(f"Unsetting `{parameter}` in {CONFIG_FILE}")
    _delete_setting(parameter, CONFIG_FILE)


def _delete_setting(key: str, config_file: str | Path) -> None:
    """
    Remove the quacc setting from the configuration file.

    Parameters
    ----------
    key
        The key in the YAML file to unset.

    Returns
    -------
    None
    """
    yaml = ruamel.yaml.YAML()
    if config_file.exists():
        with open(config_file, "r") as yaml_file:
            yaml_content = yaml.load(yaml_file)

    if yaml_content:
        yaml_content.pop(key, None)
        with open(config_file, "w") as yaml_file:
            if yaml_content:
                yaml.dump(yaml_content, yaml_file)


def _update_setting(key: str, value: Any, config_file: str | Path) -> None:
    """
    Update the quacc setting from the configuration file.

    Parameters
    ----------
    key
        The key in the YAML file to (re)set.
    value
        The value to set.

    Returns
    -------
    None
    """

    yaml = ruamel.yaml.YAML()

    if config_file.exists() and config_file.stat().st_size > 0:
        with open(config_file, "r") as yaml_file:
            yaml_content = yaml.load(yaml_file)
    else:
        yaml_content = {}

    yaml_content[key] = value
    with open(config_file, "w") as yaml_file:
        yaml.dump(yaml_content, yaml_file)


if __name__ == "__main__":
    app()
