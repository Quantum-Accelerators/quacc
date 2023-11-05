"""Atoms handling with Phonopy"""
from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure

try:
    import phonopy
except ImportError:
    phonopy = None

if TYPE_CHECKING:
    from ase import Atoms
    from numpy.typing import ArrayLike
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms


def atoms_to_phonopy(
    atoms: Atoms, supercell_matrix: ArrayLike, atom_disp: float
) -> Phonopy:
    """
    Convert an ASE atoms object to a phonopy object with displacements
    generated.

    Parameters
    ----------
    atoms
        ASE atoms object
    supercell_matrix
        Supercell matrix to use.
    atom_disp
        Atomic displacement (A).

    Returns
    -------
    Phonopy
        Phonopy object
    """
    structure = AseAtomsAdaptor().get_structure(atoms)
    phonopy_atoms = get_phonopy_structure(structure)
    phonon = phonopy.Phonopy(phonopy_atoms, supercell_matrix)
    phonon.generate_displacements(distance=atom_disp)
    return phonon


def phonopy_atoms_to_ase_atoms(phonpy_atoms: PhonopyAtoms) -> Atoms:
    """
    Convert a phonopy atoms object to an ASE atoms object.

    Parameters
    ----------
    phonpy_atoms
        Phonopy atoms object

    Returns
    -------
    Atoms
        ASE atoms object
    """
    pmg_structure = get_pmg_structure(phonpy_atoms)
    return AseAtomsAdaptor().get_atoms(pmg_structure)
