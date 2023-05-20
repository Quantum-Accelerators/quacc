"""Init data for Quacc"""
from ase import Atoms
from ase.io.jsonio import decode, encode

from quacc._version import __version__
from quacc.settings import QuaccSettings

try:
    import covalent as ct
except ImportError:
    ct = None


def atoms_as_dict(s):
    # Uses Monty's MSONable spec
    # Normally, we would want to this to be a wrapper around atoms.todict() with @module and
    # @class key-value pairs inserted. However, atoms.todict()/atoms.fromdict() does not currently
    # work properly with constraints.

    d = {"@module": "ase.atoms", "@class": "Atoms", "atoms_json": encode(s)}

    return d


def atoms_from_dict(d):
    # Uses Monty's MSONable spec
    # Normally, we would want to have this be a wrapper around atoms.fromdict()
    # that just ignores the @module/@class key-value pairs. However, atoms.todict()/atoms.fromdict()
    # does not currently work properly with constraints.

    atoms = decode(d["atoms_json"])

    return atoms


Atoms.as_dict = atoms_as_dict
Atoms.from_dict = atoms_from_dict
SETTINGS = QuaccSettings()

if ct:
    ct_config = ct.get_config()
    print("Configuring Covalent...")

    # Ensure that the create_unique_workdir is set to True
    for executor in ["dask", "local"]:
        try:
            create_unique_workdir = ct_config["executors"][executor][
                "create_unique_workdir"
            ]
        except KeyError:
            create_unique_workdir = None

        if not create_unique_workdir:
            print(f"Setting executors.{executor}.create_unique_workdir: True")
            ct.set_config({f"executors.{executor}.create_unique_workdir": True})
            ct_config = ct.get_config()

    # Ensure that use_srun is False in Slurm executor if the plugin is installed
    if "slurm" in ct_config["executors"]:
        if ct_config["executors"]["slurm"].get("use_srun", True) is not False:
            print("Setting executors.slurm.use_srun: False")
            ct.set_config({"executors.slurm.use_srun": False})
            ct_config = ct.get_config()
        if (
            ct_config["executors"]["slurm"].get("create_unique_workdir", False)
            is not True
        ):
            print("Setting executors.slurm.create_unique_workdir: True")
            ct.set_config({"executors.slurm.create_unique_workdir": True})
            ct_config = ct.get_config()

    # Make sure that the create_unique_workdir is set to True for any optional plugin
    # where this option exists, in case we missed any above
    for executor in ct_config["executors"]:
        if (
            "create_unique_workdir" in ct_config["executors"][executor]
            and ct_config["executors"][executor]["create_unique_workdir"] is not True
        ):
            print(f"Setting executors.{executor}.create_unique_workdir: True")
            ct.set_config({f"executors.{executor}.create_unique_workdir": True})
            ct_config = ct.get_config()
