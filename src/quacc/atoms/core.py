"""Utility functions for dealing with Atoms."""

from __future__ import annotations

import hashlib
import logging
from copy import deepcopy
from typing import TYPE_CHECKING

import numpy as np
from ase.filters import Filter
from ase.io.jsonio import encode
from pymatgen.io.ase import AseAtomsAdaptor

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from ase.optimize.optimize import Dynamics
    from numpy.typing import NDArray

logger = logging.getLogger(__name__)


def get_atoms_id(atoms: Atoms) -> str:
    """
    Returns a unique ID for the Atoms object. Note: The .info dict and calculator is
    excluded from the hash generation.

    Parameters
    ----------
    atoms
        Atoms object

    Returns
    -------
    str
        MD5 hash of the Atoms object
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

    return hashlib.md5(encoded_atoms.encode("utf-8"), usedforsecurity=False).hexdigest()


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


def get_charge_attribute(atoms: Atoms) -> int | None:
    """
    Get the charge of an Atoms object.

    Parameters
    ----------
    atoms
        Atoms object

    Returns
    -------
    int | None
        Charge of the Atoms object
    """
    return (
        atoms.charge  # type: ignore[attr-defined]
        if getattr(atoms, "charge", None)
        else (
            round(atoms.get_initial_charges().sum())
            if atoms.has("initial_charges")
            else None
        )
    )


def get_spin_multiplicity_attribute(atoms: Atoms) -> int | None:
    """
    Get the spin multiplicity of an Atoms object.

    Parameters
    ----------
    atoms
        Atoms object

    Returns
    -------
    int | None
        Spin multiplicity of the Atoms object
    """
    if getattr(atoms, "spin_multiplicity", None):
        return atoms.spin_multiplicity  # type: ignore[attr-defined]
    elif (
        getattr(atoms, "calc", None) is not None
        and getattr(atoms.calc, "results", None) is not None
        and atoms.calc.results.get("magmom", None) is not None
    ):
        return round(abs(atoms.calc.results["magmom"])) + 1
    elif (
        getattr(atoms, "calc", None) is not None
        and getattr(atoms.calc, "results", None) is not None
        and atoms.calc.results.get("magmoms", None) is not None
    ):
        return round(np.abs(atoms.calc.results["magmoms"].sum())) + 1
    elif atoms.has("initial_magmoms"):
        return round(np.abs(atoms.get_initial_magnetic_moments().sum())) + 1
    else:
        return None


def check_charge_and_spin(
    atoms: Atoms, charge: int | None = None, spin_multiplicity: int | None = None
) -> tuple[int, int]:
    """
    Check the validity of a given `charge` and `multiplicity`. If they are `None`, then
    set the charge and/or spin multiplicity of a molecule using the information available,
    raising a `ValueError` if there is an incompatibility.

    Parameters
    ----------
    atoms
        Atoms object
    charge
        Molecular charge
    spin_multiplicity
        Molecular multiplicity

    Returns
    -------
    charge, multiplicity
    """
    charge = charge if charge is not None else get_charge_attribute(atoms)
    spin_multiplicity = (
        spin_multiplicity
        if spin_multiplicity is not None
        else get_spin_multiplicity_attribute(atoms)
    )

    if charge is None and spin_multiplicity is not None:
        charge = 0

    try:
        mol = AseAtomsAdaptor.get_molecule(atoms)
        if charge is not None:
            if spin_multiplicity is not None:
                mol.set_charge_and_spin(charge, spin_multiplicity)
            else:
                mol.set_charge_and_spin(charge)
    except ValueError:
        mol = AseAtomsAdaptor.get_molecule(atoms, charge_spin_check=False)
        nelectrons = mol.nelectrons - charge if charge else mol.nelectrons
        default_spin_multiplicity = 1 if nelectrons % 2 == 0 else 2
        mol.set_charge_and_spin(
            charge if charge is not None else mol.charge,
            (
                spin_multiplicity
                if spin_multiplicity is not None
                else default_spin_multiplicity
            ),
        )
    if (mol.nelectrons + mol.spin_multiplicity) % 2 != 1:
        raise ValueError(
            f"Charge of {mol.charge} and spin multiplicity of {mol.spin_multiplicity} is"
            " not possible for this molecule."
        )
    logger.debug(
        f"Setting charge to {mol.charge} and spin multiplicity to {mol.spin_multiplicity}"
    )

    return mol.charge, mol.spin_multiplicity


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
        Nx3 matrix, where N is the number of atoms. This means that there is potentially a different translation
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
