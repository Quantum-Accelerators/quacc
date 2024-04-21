"""Atoms handling with Phonopy."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
from monty.dev import requires
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

try:
    import phonopy

    has_deps = find_spec("seekpath") is not None
except ImportError:
    has_deps = False

if TYPE_CHECKING:
    from ase.atoms import Atoms

    if phonopy:
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms


@requires(has_deps, "Phonopy or seekpath is not installed.")
def prep_phonopy(
    atoms: Atoms,
    fixed_atoms: list[int] | None = None,
    min_lengths: float | tuple[float, float, float] | None = None,
    phonopy_kwargs: dict | None = None,
    generate_displacements_kwargs: dict | None = None,
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
    generate_displacements_kwargs = generate_displacements_kwargs or {}

    symprec = phonopy_kwargs.pop("symprec", 1e-5)

    supercell_matrix = phonopy_kwargs.pop("supercell_matrix", None)

    structure = AseAtomsAdaptor.get_structure(atoms)
    structure = SpacegroupAnalyzer(
        structure, symprec=symprec
    ).get_symmetrized_structure()
    atoms = AseAtomsAdaptor.get_atoms(structure)

    fixed, non_fixed = atoms[fixed_atoms], atoms[~fixed_atoms]

    non_fixed = AseAtomsAdaptor.get_structure(non_fixed)

    if supercell_matrix is None and min_lengths is not None:
        n_supercells = np.round(np.ceil(min_lengths / atoms.cell.lengths()))
        supercell_matrix = np.diag([n_supercells, n_supercells, n_supercells])

    phonopy_atoms = get_phonopy_structure(non_fixed)
    phonon = phonopy.Phonopy(
        phonopy_atoms,
        symprec=symprec,
        supercell_matrix=supercell_matrix,
        **phonopy_kwargs,
    )
    phonon.generate_displacements(**generate_displacements_kwargs)

    if len(fixed) > 0:
        fixed = get_phonopy_structure(AseAtomsAdaptor.get_structure(non_fixed))
        fixed = phonopy.structure.cells.get_supercell(non_fixed, supercell_matrix)
        fixed = AseAtomsAdaptor.get_atoms(get_pmg_structure(non_fixed))

    return phonon, fixed


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
