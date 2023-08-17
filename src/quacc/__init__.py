"""Init data for quacc"""

import os
from importlib.metadata import version

from ase import Atoms
from ase.io.jsonio import decode, encode

from quacc.settings import QuaccSettings
from quacc.util.wflows import flow, job, subflow


def atoms_as_dict(s):
    # Uses Monty's MSONable spec
    # Normally, we would want to this to be a wrapper around atoms.todict() with @module and
    # @class key-value pairs inserted. However, atoms.todict()/atoms.fromdict() does not currently
    # work properly with constraints.
    return {"@module": "ase.atoms", "@class": "Atoms", "atoms_json": encode(s)}


def atoms_from_dict(d):
    # Uses Monty's MSONable spec
    # Normally, we would want to have this be a wrapper around atoms.fromdict()
    # that just ignores the @module/@class key-value pairs. However, atoms.todict()/atoms.fromdict()
    # does not currently work properly with constraints.
    return decode(d["atoms_json"])


# Load the version
__version__ = version("quacc")

# Make Atoms MSONable
Atoms.as_dict = atoms_as_dict
Atoms.from_dict = atoms_from_dict

# Load the settings
SETTINGS = QuaccSettings()

for f in {SETTINGS.SCRATCH_DIR, SETTINGS.RESULTS_DIR}:
    if f:
        os.makedirs(f, exist_ok=True)
