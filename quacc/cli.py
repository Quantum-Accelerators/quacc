"""Quacc command-line interface"""
import sys


def main() -> None:
    """
    Set up command-line interface.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """
    if len(sys.argv) > 1 and sys.argv[1] == "config":
        configure_covalent()
    else:
        print("Invalid command.")


def configure_covalent() -> None:
    """
    Configure Covalent to work with Quacc.

    This function is called by the `quacc config` command.

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


if __name__ == "__main__":
    main()
