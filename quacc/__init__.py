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
    if "slurm" in ct_config["executors"]:
        if ct_config["executors"]["slurm"].get("use_srun", True) is not False:
            print("Configuring Covalent... Setting executors.slurm.use_srun: False")
            ct.set_config({"executors.slurm.use_srun": False})
            ct_config = ct.get_config()
