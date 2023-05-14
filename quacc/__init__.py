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
    for executor in ct_config["executors"]:
        if ct_config["executors"][executor].get("create_unique_workdir", "false") == "false":
            warnings.warn(
                f"Updating Covalent configuration... setting executors.{executor}.create_unique_workdir: 'true'",
                UserWarning,
            )
            ct.set_config({f"executors.{executor}.create_unique_workdir": "true"})
    if ct_config["executors"].get("slurm", None):
        if ct_config["executors"]["slurm"].get("use_srun", "true") == "true":
            warnings.warn(
                "Updating Covalent configuration... setting executors.slurm.use_srun: 'false'",
                UserWarning,
            )
            ct.set_config({"executors.slurm.use_srun": "false"})
