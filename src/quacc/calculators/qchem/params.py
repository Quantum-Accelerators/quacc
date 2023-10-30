"""
Parameter-related utilities for the Q-Chem calculator.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor

if TYPE_CHECKING:
    from typing import Any, Literal

    from pymatgen.core.structure import Molecule

logger = logging.getLogger(__name__)


def get_rem_swaps(rem: dict[str, Any]) -> dict[str, Any]:
    """
    Automatic swaps for the rem dictionary.

    Parameters
    ----------
    rem
        rem dictionary

    Returns
    -------
    dict
        rem dictionary with swaps
    """
    if "scf_guess" not in rem:
        logger.info("Copilot: Setting scf_guess in `rem` to 'read'")
        rem["scf_guess"] = "read"
    if "max_scf_cycles" not in rem and rem.get("scf_algorithm") == "gdm":
        logger.info("Copilot: Setting max_scf_cycles in `rem` to 200")
        rem["max_scf_cycles"] = 200

    return rem


def get_molecule(
    atoms: Atoms | list[Atoms] | Literal["read"],
    charge: int,
    spin_multiplicity: int,
) -> Molecule | list[Molecule] | Literal["read"]:
    """
    Convert ASE Atom(s) to Molecule(s) suitable for `QCInput`.

    Parameters
    ----------
    atoms
        Input `atoms` kwarg to the calculator.

    Returns
    -------
    Molecule | list[Molecule] | Literal["read"]
        The corresponding `molecule` kwarg to pass to QCInput.
    """
    adaptor = AseAtomsAdaptor()

    if isinstance(atoms, Atoms):
        pmg_obj = adaptor.get_molecule(atoms)
        pmg_obj.set_charge_and_spin(charge, spin_multiplicity)
        return pmg_obj

    if isinstance(atoms, list):
        molecules = []
        for atoms_ in atoms:
            pmg_obj = adaptor.get_molecule(atoms_)
            pmg_obj.set_charge_and_spin(charge, spin_multiplicity)
            molecules.append(pmg_obj)
        return molecules

    if isinstance(atoms, str):
        return atoms

    raise TypeError(f"Invalid type for atoms: {type(atoms)}")
