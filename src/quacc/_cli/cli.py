"""Quacc CLI module."""
from __future__ import annotations

import click


@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx: click.Context) -> None:
    """
    Welcome to the quacc command-line interface.

    Parameters
    ----------
    ctx
        Click context object

    Returns
    -------
    None
    """

    # Return help message if no command is provided
    if ctx.invoked_subcommand is None:
        ctx = click.get_current_context()
        click.echo(ctx.get_help())


if __name__ == "__main__":
    cli()
