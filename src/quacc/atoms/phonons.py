"""Atoms handling with Phonopy."""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from monty.dev import requires
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

try:
    import phonopy
except ImportError:
    phonopy = None

if TYPE_CHECKING:
    from ase.atoms import Atoms

    if phonopy:
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms


@requires(phonopy, "Phonopy is not installed.")
def get_phonopy(
    atoms: Atoms,
    min_length: float | None = None,
    symprec: float = 1e-5,
    displacement: float = 0.01,
    phonopy_kwargs: dict | None = None,
) -> Phonopy:
    """
    Convert an ASE atoms object to a phonopy object with displacements generated.

    Parameters
    ----------
    atoms
        ASE atoms object.
    min_length
        Minimum length of each lattice dimension (A).
    symprec
        Precision for symmetry detection.
    displacement
        Atomic displacement (A).
    phonopy_kwargs
        Additional kwargs to pass to the Phonopy class.

    Returns
    -------
    Phonopy
        Phonopy object
    """
    phonopy_kwargs = phonopy_kwargs or {}

    structure = AseAtomsAdaptor().get_structure(atoms)
    structure = SpacegroupAnalyzer(
        structure, symprec=symprec
    ).get_symmetrized_structure()

    if min_length:
        n_supercells = np.round(np.ceil(min_length / atoms.cell.lengths()))
        supercell_matrix = np.diag([n_supercells, n_supercells, n_supercells])
    else:
        supercell_matrix = None

    phonopy_atoms = get_phonopy_structure(structure)
    phonon = phonopy.Phonopy(
        phonopy_atoms,
        symprec=symprec,
        supercell_matrix=supercell_matrix,
        **phonopy_kwargs,
    )
    phonon.generate_displacements(distance=displacement)
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
