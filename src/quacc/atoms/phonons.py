"""Atoms handling with Phonopy."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
from monty.dev import requires
from pymatgen.core import Structure
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure

has_phonopy = bool(find_spec("phonopy"))

if has_phonopy:
    from phonopy import Phonopy

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from numpy.typing import NDArray

    if has_phonopy:
        from phonopy.structure.atoms import PhonopyAtoms


@requires(has_phonopy, "Phonopy not installed.")
def get_phonopy(
    atoms: Atoms,
    min_lengths: float | tuple[float, float, float] | None = None,
    supercell_matrix: (
        NDArray
        | tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
        | None
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

    if supercell_matrix is None and min_lengths is not None:
        supercell_matrix = np.diag(
            np.round(np.ceil(min_lengths / atoms.cell.lengths()))
        )

    phonon = Phonopy(
        get_phonopy_structure(Structure.from_ase_atoms(atoms)),
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
    return pmg_structure.to_ase_atoms()
