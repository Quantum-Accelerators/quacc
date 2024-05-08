"""Atoms handling with Phonopy."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ase.constraints import FixAtoms
from monty.dev import requires
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

try:
    import phonopy

    has_phonopy = True
except ImportError:
    has_phonopy = False

if TYPE_CHECKING:
    from ase.atoms import Atoms

    if phonopy:
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms


@requires(has_phonopy, "Phonopy not installed.")
def get_phonopy(
    atoms: Atoms,
    min_lengths: float | tuple[float, float, float] | None = None,
    supercell_matrix: (
        tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]] | None
    ) = None,
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
    min_lengths
        Minimum length of each lattice dimension (A).
    supercell_matrix
        The supercell matrix to use. If specified, it will override any
        value specified by `min_lengths`.
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

    fixed_indices = np.array(
        [
            constr.get_indices()
            for constr in atoms.constraints
            if isinstance(constr, FixAtoms)
        ]
    )
    
    fixed_indices.flatten()
    is_fixed_atoms = np.array([i in fixed_indices for i in range(len(atoms))])

    structure = AseAtomsAdaptor.get_structure(atoms)
    structure = SpacegroupAnalyzer(
        structure, symprec=symprec
    ).get_symmetrized_structure()
    atoms = structure.to_ase_atoms()

    fixed_atoms, non_fixed_atoms = atoms[is_fixed_atoms], atoms[~is_fixed_atoms]

    non_fixed_atoms = AseAtomsAdaptor.get_structure(non_fixed_atoms)

    if supercell_matrix is None and min_lengths is not None:
        supercell_matrix = np.diag(np.round(np.ceil(min_lengths / atoms.cell.lengths())))

    phonopy_atoms = get_phonopy_structure(non_fixed_atoms)
    phonon = phonopy.Phonopy(
        phonopy_atoms,
        symprec=symprec,
        supercell_matrix=supercell_matrix,
        **phonopy_kwargs,
    )
    phonon.generate_displacements(distance=displacement)

    if fixed_atoms:
        fixed_atoms = phonopy_atoms_to_ase_atoms(
            phonopy.structure.cells.get_supercell(
                get_phonopy_structure(AseAtomsAdaptor.get_structure(fixed_atoms)),
                supercell_matrix,
            )
        )

    return phonon, fixed_atoms


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
    return pmg_structure.to_ase_atoms()
