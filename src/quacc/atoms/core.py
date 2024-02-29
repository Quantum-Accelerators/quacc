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

    return hashlib.md5(encoded_atoms.encode("utf-8")).hexdigest()


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
    if atoms.pbc.any():
        struct = AseAtomsAdaptor.get_structure(atoms)
    else:
        struct = AseAtomsAdaptor.get_molecule(atoms, charge_spin_check=False)

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
        atoms.charge
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
    return (
        atoms.spin_multiplicity
        if getattr(atoms, "spin_multiplicity", None)
        else (
            round(np.abs(atoms.get_initial_magnetic_moments().sum()) + 1)
            if atoms.has("initial_magmoms")
            else None
        )
    )


def check_charge_and_spin(
    atoms: Atoms, charge: int | None = None, spin_multiplicity: int | None = None
) -> tuple[int, int]:
    """
    Check the validity of a given `charge` and `multiplicity`. If they are `None`, then
    set the charge and/or spin multiplicity of a molecule using the following routine,
    raising a `ValueError` if there is an incompatibility.

    Charges:

    1. If `charge` is specified, that is the charge.

    2. If `atoms.charge` is present, that is the charge.

    3. If `atoms.has("initial_charges")`, then
    `atoms.get_initial_charges.sum()` is the charge.

    4. If `spin_multiplicity` is set, set the charge to the smallest physically
    possible value.

    5. Otherwise, set to 0.

    Spin multiplicity:

    1. If `spin_multiplicity` is specified, that is the spin multiplicity.

    2. If `atoms.spin_multiplicity` is present, that is the spin multiplicity.

    3. If `atoms.has("initial_magmoms")`, then
    `np.abs(atoms.get_initial_magnetic_moments().sum())+1` is the spin
    multiplicity.

    4. If none of the above, use Pymatgen to identify the lowest physically
    possible spin multiplicity given the number of electrons and the charge, if
    set.

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


def get_final_atoms_from_dyn(dyn: Dynamics) -> Atoms:
    """
    Get the final atoms object from a dynamics run.

    Parameters
    ----------
    dyn
        ASE dynamics object

    Returns
    -------
    Atoms
        Atoms object
    """
    return dyn.atoms.atoms if isinstance(dyn.atoms, Filter) else dyn.atoms
