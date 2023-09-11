"""
Utility functions for dealing with Atoms
"""
from __future__ import annotations

import hashlib
import logging
from copy import deepcopy
from typing import TYPE_CHECKING

import numpy as np
from ase.io.jsonio import encode
from pymatgen.io.ase import AseAtomsAdaptor

if TYPE_CHECKING:
    from ase import Atoms

logger = logging.getLogger(__name__)


def prep_next_run(
    atoms: Atoms, assign_id: bool = True, move_magmoms: bool = True
) -> Atoms:
    """
    Prepares the Atoms object for a new run.

    Depending on the arguments, this function will:
        - Move the converged magnetic moments to the initial magnetic moments.
        - Assign a unique ID to the Atoms object in atoms.info["_id"]. Any
          existing IDs will be moved to atoms.info["_old_ids"].

    In all cases, the calculator will be reset so new jobs can be run.

    Parameters
    ----------
    atoms
        Atoms object
    assign_id
        Whether to assign a unique ID to the Atoms object in atoms.info["_id"].
        Any existing IDs will be moved to atoms.info["_old_ids"].
    move_magmoms
        If True, move atoms.calc.results["magmoms"] to
        atoms.get_initial_magnetic_moments()

    Returns
    -------
    Atoms
        Updated Atoms object.
    """
    atoms = copy_atoms(atoms)

    if (
        move_magmoms
        and hasattr(atoms, "calc")
        and getattr(atoms.calc, "results", None) is not None
    ):
        # If there are initial magmoms set, then we should see what the final
        # magmoms are. If they are present, move them to initial. If they are
        # not present, it means the calculator doesn't support the "magmoms"
        # property so we have to retain the initial magmoms given no further
        # info.
        if atoms.has("initial_magmoms"):
            atoms.set_initial_magnetic_moments(
                atoms.calc.results.get("magmoms", atoms.get_initial_magnetic_moments())
            )
        # If there are no initial magmoms set, just check the results and set
        # everything to 0.0 if there is nothing there.
        else:
            atoms.set_initial_magnetic_moments(
                atoms.calc.results.get("magmoms", [0.0] * len(atoms))
            )

    # Clear off the calculator so we can run a new job. If we don't do this,
    # then something like atoms *= (2,2,2) still has a calculator attached,
    # which is a bit confusing.
    atoms.calc = None

    # Give the Atoms object a unique ID. This will be helpful for querying
    # later. Also store any old IDs somewhere else for future reference. Note:
    # Keep this at the end of the function so that the ID is assigned based on
    # the returned Atoms object.
    if assign_id:
        if atoms.info.get("_id", None) is not None:
            if atoms.info.get("_old_ids") is None:
                atoms.info["_old_ids"] = []
            atoms.info["_old_ids"].append(atoms.info["_id"])
        atoms.info["_id"] = get_atoms_id(atoms)

    return atoms


def set_magmoms(
    atoms: Atoms,
    elemental_mags_dict: dict | None = None,
    elemental_mags_default: float = 1.0,
    copy_magmoms: bool = True,
    mag_cutoff: float | None = 0.05,
) -> Atoms:  # sourcery skip
    """
    Sets the initial magnetic moments in the Atoms object.

    This function deserves particular attention. The following logic is applied:
    - If there is a converged set of magnetic moments, those are moved to the
    initial magmoms if copy_magmoms is True. - If there is no converged set of
    magnetic moments but the user has set initial magmoms, those are simply used
    as is. - If there are no converged magnetic moments or initial magnetic
    moments, then the default magnetic moments from the preset
    elemental_mags_dict (if specified) are set as the initial magnetic moments.
    - For any of the above scenarios, if mag_cutoff is not None, the newly set
    initial magnetic moments are checked. If all have a magnitude below
    mag_cutoff, then they are all set to 0 (no spin polarization).

    Parameters
    ----------
    atoms
        Atoms object
    elemental_mags_dict
        Dictionary of elemental symbols and their corresponding magnetic moments
        to set. If None, no default values will be used.
    elemental_mags_default
        Default magnetic moment on an element if no magnetic moment is specified
        in the elemental_mags_dict. Only used if elemental_mags_dict is not
        None. This kwarg is mainly a convenience so that you don't need to list
        every single element in the elemental_mags_dict.
    copy_magmoms
        Whether to copy the magnetic moments from the converged set of magnetic
        moments to the initial magnetic moments.
    mag_cutoff
        Magnitude below which the magnetic moments are considered to be zero. If
        None, no cutoff will be applied

    Returns
    -------
    Atoms
        Atoms object
    """

    # Handle the magnetic moments Check if a prior job was run and pull the
    # prior magmoms
    if hasattr(atoms, "calc") and getattr(atoms.calc, "results", None) is not None:
        mags = atoms.calc.results.get("magmoms", [0.0] * len(atoms))
        # Note: It is important that we set mags to 0.0 here rather than None if
        # the calculator has no magmoms because: 1) ispin=1 might be set, and 2)
        # we do not want the preset magmoms to be used.
    else:
        mags = None

    # Check if the user has set any initial magmoms
    has_initial_mags = atoms.has("initial_magmoms")

    # If there are no initial magmoms set and this is not a follow-up job, we
    # may need to add some from the preset yaml.
    if mags is None:
        if not has_initial_mags:
            # If the preset dictionary has default magmoms, set those by
            # element. If the element isn't in the magmoms dict then set it to
            # mag_default.
            if elemental_mags_dict:
                initial_mags = np.array(
                    [
                        elemental_mags_dict.get(atom.symbol, elemental_mags_default)
                        for atom in atoms
                    ]
                )
                atoms.set_initial_magnetic_moments(initial_mags)
        else:
            pass
    elif copy_magmoms:
        atoms.set_initial_magnetic_moments(mags)

    # If all the set mags are below mag_cutoff, set them to 0
    if mag_cutoff:
        has_new_initial_mags = atoms.has("initial_magmoms")
        new_initial_mags = atoms.get_initial_magnetic_moments()
        if has_new_initial_mags and np.all(np.abs(new_initial_mags) < mag_cutoff):
            atoms.set_initial_magnetic_moments([0.0] * len(atoms))

    return atoms


def get_atoms_id(atoms: Atoms) -> str:
    """
    Returns a unique ID for the Atoms object. Note: The .info dict and
    calculator is excluded from the hash generation.

    Parameters
    ----------
    atoms
        Atoms object

    Returns
    -------
    md5hash
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
        struct = AseAtomsAdaptor.get_molecule(atoms)

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


def check_charge_and_spin(
    atoms: Atoms,
    charge: int | None = None,
    spin_multiplicity: int | None = None,
) -> (int, int):
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

    charge = (
        charge
        if charge is not None
        else atoms.charge
        if getattr(atoms, "charge", None)
        else round(atoms.get_initial_charges().sum())
        if atoms.has("initial_charges")
        else None
    )

    spin_multiplicity = (
        spin_multiplicity
        if spin_multiplicity is not None
        else atoms.spin_multiplicity
        if getattr(atoms, "spin_multiplicity", None)
        else round(np.abs(atoms.get_initial_magnetic_moments().sum()) + 1)
        if atoms.has("initial_magmoms")
        else None
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
            spin_multiplicity
            if spin_multiplicity is not None
            else default_spin_multiplicity,
        )
    if (mol.nelectrons + mol.spin_multiplicity) % 2 != 1:
        raise ValueError(
            f"Charge of {mol.charge} and spin multiplicity of {mol.spin_multiplicity} is"
            " not possible for this molecule."
        )
    logger.info(
        f"Setting charge to {mol.charge} and spin multiplicity to {mol.spin_multiplicity}"
    )

    return mol.charge, mol.spin_multiplicity
