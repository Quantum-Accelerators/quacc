from ase.atoms import Atom, Atoms
from ase.io.jsonio import decode, encode

from quacc._version import __version__


def atoms_as_dict(self):
    # Uses Monty's MSONable spec
    # Normally, we would want to this to be a wrapper around atoms.todict() with @module and
    # @class key-value pairs inserted. However, atoms.todict()/atoms.fromdict() does not currently
    # work properly with constraints.

    d = {"@module": "ase.atoms", "@class": "Atoms", "atoms_json": encode(self)}

    return d


def atoms_from_dict(d):
    # Uses Monty's MSONable spec
    # Normally, we would want to have this be a wrapper around atoms.fromdict()
    # that just ignores the @module/@class key-value pairs. However, atoms.todict()/atoms.fromdict()
    # does not currently work properly with constraints.

    atoms = decode(d["atoms_json"])

    return atoms


def atom_as_dict(self):
    # Uses Monty's MSONable spec
    d = {"@module": "ase.atoms", "@class": "Atom", "atom_json": encode(self)}

    return d


def atom_from_dict(self, d):
    # Uses Monty's MSONable spec
    atom = decode(d["atom_json"])

    return atom


Atoms.as_dict = atoms_as_dict
Atoms.from_dict = atoms_from_dict
Atom.as_dict = atom_as_dict
Atom.from_dict = atom_from_dict
