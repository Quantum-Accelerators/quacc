"""
Parameter-related utilities for the Q-Chem calculator.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.atoms.core import atoms_to_pmg

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase import Atoms
    from pymatgen.core.structure import Molecule


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
        rem["scf_guess"] = "read"
    if "max_scf_cycles" not in rem and rem.get("scf_algorithm") == "gdm":
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
    molecule
        The corresponding `molecule` kwarg to pass to QCInput.
    """

    if isinstance(atoms, Atoms):
        return atoms_to_pmg(atoms, charge=charge, spin_multiplicity=spin_multiplicity)
    if isinstance(atoms, list):
        return [
            atoms_to_pmg(
                atoms_, charge=charge, spin_multiplicity=spin_multiplicity
            )
            for atoms_ in atoms
        ]
    if isinstance(atoms, str):
        return atoms
