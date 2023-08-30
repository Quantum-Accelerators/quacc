"""Quacc CLI module."""
from __future__ import annotations

import os
from typing import TYPE_CHECKING

import ruamel.yaml
import typer

from quacc import SETTINGS

if TYPE_CHECKING:
    from typing import Any

app = typer.Typer()


CONFIG_FILE = SETTINGS.CONFIG_FILE or os.path.expanduser("~/.quacc.yaml")


@app.command()
def set(parameter: str, new_value: str) -> None:
    """
    Set the quacc variable.

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
    typer.echo(f"Setting `{parameter}` to `{new_value}` in {CONFIG_FILE}")
    _update_settings(key=parameter, value=new_value)


def _update_settings(key: str | None = None, value: Any | None = None) -> None:
    """
    Update the quacc settings file.

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

    if not os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, "w") as f:
            f.write("")

    # Load the YAML content from the file
    yaml = ruamel.yaml.YAML()
    with open(CONFIG_FILE, "r") as yaml_file:
        yaml_content = yaml.load(yaml_file)

    if yaml_content and key in yaml_content:
        yaml_content[key] = value
        del yaml_content[key]

    # Write the updated content back to the file
    with open(CONFIG_FILE, "w") as yaml_file:
        yaml.dump(yaml_content, yaml_file)
