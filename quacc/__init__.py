"""Init data for QuAcc"""
import warnings

from ase.atoms import Atoms
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
    # Ensure that the create_unique_workdir is set to True
    # These are the executors where we know this parameter exists
    for executor in ["dask", "local", "slurm"]:
        try:
            create_unique_workdir = ct.get_config()["executors"][executor][
                "create_unique_workdir"
            ]
        except:
            create_unique_workdir = None

        if not create_unique_workdir:
            warnings.warn(
                f"Updating Covalent configuration... setting executors.{executor}.create_unique_workdir: True",
                UserWarning,
            )
            ct.set_config({f"executors.{executor}.create_unique_workdir": True})

    # Make sure that the create_unique_workdir is set to True for any executor
    # where this option exists, in case we missed any above
    for executor in ct.get_config()["executors"]:
        if (
            "create_unique_workdir" in ct.get_config()["executors"][executor]
            and ct.get_config()["executors"][executor]["create_unique_workdir"]
            is not True
        ):
            warnings.warn(
                f"Updating Covalent configuration... setting executors.{executor}.create_unique_workdir: True",
                UserWarning,
            )
            ct.set_config({f"executors.{executor}.create_unique_workdir": True})

    # Ensure that use_srun is False in Slurm executor if the plugin is installed
    if "slurm" in ct.get_config()["executors"]:
        if ct.get_config()["executors"]["slurm"].get("use_srun", True) is not False:
            warnings.warn(
                "Updating Covalent configuration... setting executors.slurm.use_srun: False",
                UserWarning,
            )
            ct.set_config({"executors.slurm.use_srun": False})
