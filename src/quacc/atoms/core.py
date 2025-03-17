"""Utility functions for dealing with Atoms."""

from __future__ import annotations

from copy import deepcopy
from hashlib import md5
from logging import getLogger
from typing import TYPE_CHECKING

import numpy as np
from ase.filters import Filter
from ase.io.jsonio import encode
from pymatgen.io.ase import AseAtomsAdaptor

if TYPE_CHECKING:
    from hashlib import _Hash

    from ase.atoms import Atoms
    from ase.optimize.optimize import Dynamics
    from numpy.typing import NDArray

LOGGER = getLogger(__name__)


def _encode_atoms(atoms: Atoms) -> _Hash:
    """
    Returns a byte encoding for the Atoms object. Note: The .info dict and calculator is excluded.

    Parameters
    ----------
    atoms
        Atoms object

    Returns
    -------
    _Hash
        Encoded Atoms object
    """
    atoms = copy_atoms(atoms)
    atoms.info = {}
    atoms.calc = None
    encoded_atoms = encode(atoms)
    # This is a hack to avoid int32/int64 and float32/float64 differences
    # between machines.
    encoded_atoms = (
        encoded_atoms.replace("int64", "int")
        .replace("int32", "int")
        .replace("float64", "float")
        .replace("float32", "float")
    )

    return md5(encoded_atoms.encode("utf-8"), usedforsecurity=False)


def get_atoms_id(atoms: Atoms) -> str:
    """
    Get a unique identifier for an Atoms object.

    Parameters
    ----------
    atoms
        Atoms object

    Returns
    -------
    str
        Unique identifier for the Atoms object in the form of a string
    """
    return _encode_atoms(atoms).hexdigest()


def get_atoms_id_parsl(atoms: Atoms, output_ref: bool = False) -> bytes:  # noqa: ARG001
    """
    Get a Parsl compatible unique identifier for an Atoms object.

    Parameters
    ----------
    atoms
        Atoms object
    output_ref
        Parsl specific parameter, needed for Parsl to work, unused.

    Returns
    -------
    bytes
        Unique identifier for the Atoms object in the form of bytes
    """
    return _encode_atoms(atoms).digest()


def check_is_metal(atoms: Atoms) -> bool:
    """
    Checks if a structure is a likely metal.

    Parameters
    ----------
    atoms
        Atoms object

    Returns
    -------
    bool
        True if the structure is likely a metal; False otherwise
    """
    struct = (
        AseAtomsAdaptor().get_structure(atoms)
        if atoms.pbc.any()
        else AseAtomsAdaptor().get_molecule(atoms, charge_spin_check=False)
    )

    return all(k.is_metal for k in struct.composition)


def copy_atoms(atoms: Atoms) -> Atoms:
    """
    Simple function to copy an atoms object to prevent mutability.

    Parameters
    ----------
    atoms
        Atoms object

    Returns
    -------
    atoms
        Atoms object
    """
    try:
        atoms = deepcopy(atoms)
    except Exception:
        # Needed because of ASE issue #1084
        calc = atoms.calc
        atoms = atoms.copy()
        atoms.calc = calc

    return atoms


def get_spin_multiplicity_attribute(atoms: Atoms) -> int | None:
    """
    Get the spin multiplicity of an Atoms object.

    Parameters
    ----------
    atoms
        Atoms object

    Returns
    -------
    int
        Spin multiplicity of the Atoms object
    """
    if getattr(atoms, "spin_multiplicity", None):
        return atoms.spin_multiplicity  # type: ignore[attr-defined]

    try:
        results = atoms.calc.results  # type: ignore[attr-defined]
    except AttributeError:
        results = None
    if results:
        if results.get("magmom", None) is not None:
            return round(abs(results["magmom"])) + 1
        if results.get("magmoms", None) is not None:
            return round(np.abs(results["magmoms"].sum())) + 1

    if atoms.has("initial_magmoms"):
        return round(np.abs(atoms.get_initial_magnetic_moments().sum())) + 1

    LOGGER.warning("Could not determine spin multiplicity. Assuming 1.")
    return 1


def get_final_atoms_from_dynamics(dynamics: Dynamics | Filter) -> Atoms:
    """
    Get the final atoms object from a dynamics run.

    Parameters
    ----------
    dynamics
        ASE dynamics object

    Returns
    -------
    Atoms
        Atoms object
    """
    return (
        dynamics.atoms.atoms if isinstance(dynamics.atoms, Filter) else dynamics.atoms
    )


def perturb(mol: Atoms, matrix: list[list[float]] | NDArray, scale: float) -> Atoms:
    """
    Perturb each atom in a molecule by a (scaled) 1x3 vector, reflecting e.g. a vibrational normal mode.

    Parameters
    ----------
    mol
        ASE Atoms object representing a molecule
    matrix
        N x 3 matrix, where N is the number of atoms. This means that there is potentially a different translation
        vector for each atom in the molecule.
    scale
        Scaling factor for perturbation

    Returns
    -------
    Atoms
        The input molecule after perturbation
    """

    mol_copy = copy_atoms(mol)
    mode = np.asarray(matrix)

    orig_pos = mol_copy.get_positions()

    pos = orig_pos + mode * scale
    mol_copy.set_positions(pos)

    return mol_copy
