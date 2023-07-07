"""Quacc CLI module."""

import click


@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx: click.Context) -> None:
    """Welcome to the quacc command-line interface."""

    # Return help message if no command is provided
    if ctx.invoked_subcommand is None:
        ctx = click.get_current_context()
        click.echo(ctx.get_help())


@click.command(name="config")
def config_command() -> None:
    """
    Modify the Covalent configuration to work with quacc.
    Parameters
    ----------
    None
    Returns
    -------
    None
    """
    import covalent as ct

    ct_config = ct.get_config()

    # Make sure that the create_unique_workdir is set to True for any plugin
    # where this option exists
    for executor in ct_config["executors"]:
        if (
            "create_unique_workdir" in ct_config["executors"][executor]
            and ct_config["executors"][executor]["create_unique_workdir"] is not True
        ):
            print(
                f"Configuring Covalent... Setting executors.{executor}.create_unique_workdir: True"
            )
            ct.set_config({f"executors.{executor}.create_unique_workdir": True})
            ct_config = ct.get_config()

    # Ensure that use_srun is False in Slurm executor if the plugin is installed
    if (
        "slurm" in ct_config["executors"]
        and ct_config["executors"]["slurm"].get("use_srun", True) is not False
    ):
        print("Configuring Covalent... Setting executors.slurm.use_srun: False")
        ct.set_config({"executors.slurm.use_srun": False})
        ct_config = ct.get_config()

    print(
        "Covalent configuration: Complete!\nIf you clear your Covalent configuration, make sure to re-run `quacc config`."
    )


cli.add_command(config_command)

if __name__ == "__main__":
    cli()
