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
    ct_config = ct.get_config()
    if not ct_config["executors"]["local"].get("create_unique_workdir", False):
        warnings.warn(
            "Updating Covalent configuration: setting covalent.executors.local.create_unique_workdir to True",
            UserWarning,
        )
        ct.set_config({"executors.local.create_unique_workdir": True})
    if not ct_config["executors"]["dask"].get("create_unique_workdir", False):
        warnings.warn(
            "Updating Covalent configuration: setting covalent.executors.dask.create_unique_workdir to True",
            UserWarning,
        )
        ct.set_config({"executors.dask.create_unique_workdir": True})
    if ct_config["executors"].get("slurm", None):
        if not ct_config["executors"]["slurm"].get("create_unique_workdir", False):
            warnings.warn(
                "Updating Covalent configuration: setting covalent.executors.slurm.create_unique_workdir to True",
                UserWarning,
            )
            ct.set_config({"executors.slurm.create_unique_workdir": True})
        if ct_config["executors"]["slurm"].get("use_srun", True):
            warnings.warn(
                "Updating Covalent configuration: setting covalent.executors.slurm.use_srun to False",
                UserWarning,
            )
            ct.set_config({"executors.slurm.use_srun": False})
