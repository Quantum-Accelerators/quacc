"""Atoms handling with Phonopy."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
from monty.dev import requires
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

has_phonopy = find_spec("phonopy")

if has_phonopy:
    from phonopy import Phonopy

if TYPE_CHECKING:
    from ase.atoms import Atoms

    if has_phonopy:
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
    fixed_indices: list[int] | None = None,
    phonopy_kwargs: dict | None = None,
) -> tuple[Phonopy, Phonopy | None]:
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
    fixed_indices
        Indices of atoms to fix for the phonon calculation.
    phonopy_kwargs
        Additional kwargs to pass to the Phonopy class.

    Returns
    -------
    Phonopy
        Phonopy object with displacements generated for the unfixed atoms.
    Phonopy | None
        Phonopy object for the fixed atoms if `fixed_indices` is specified.
        Otherwise, None.
    """
    phonopy_kwargs = phonopy_kwargs or {}
    fixed_indices = fixed_indices or []

    symmetrized_structure = SpacegroupAnalyzer(
        AseAtomsAdaptor().get_structure(atoms), symprec=symprec
    ).get_symmetrized_structure()

    if supercell_matrix is None and min_lengths is not None:
        supercell_matrix = np.diag(
            np.round(np.ceil(min_lengths / np.array(symmetrized_structure.lattice.abc)))
        )

    unfixed_structure = symmetrized_structure.copy().remove_sites(fixed_indices)
    fixed_structure = symmetrized_structure.copy().remove_sites(
        [i for i in range(len(symmetrized_structure)) if i not in fixed_indices]
    )

    unfixed_phonopy = Phonopy(
        get_phonopy_structure(unfixed_structure),
        symprec=symprec,
        supercell_matrix=supercell_matrix,
        **phonopy_kwargs,
    )
    unfixed_phonopy.generate_displacements(distance=displacement)
    fixed_phonopy = (
        Phonopy(get_phonopy_structure(fixed_structure), supercell_matrix)
        if fixed_indices
        else None
    )

    return unfixed_phonopy, fixed_phonopy


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
