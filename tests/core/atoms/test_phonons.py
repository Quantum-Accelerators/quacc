from __future__ import annotations

import pytest

pytest.importorskip("phonopy")
pytest.importorskip("seekpath")

import numpy as np
from ase.build import bulk
from numpy.testing import assert_almost_equal, assert_array_equal

from quacc.atoms.phonons import get_atoms_supercell_by_phonopy, get_phonopy


def test_get_phonopy():
    atoms = bulk("Cu")
    phonopy = get_phonopy(atoms)
    assert_array_equal(phonopy.supercell_matrix, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    phonopy = get_phonopy(atoms, min_lengths=5)
    assert_array_equal(phonopy.supercell_matrix, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])

    phonopy = get_phonopy(atoms, min_lengths=[5, 10, 5])
    assert_array_equal(phonopy.supercell_matrix, [[2, 0, 0], [0, 4, 0], [0, 0, 2]])

    phonopy = get_phonopy(atoms, displacement=1)
    assert_almost_equal(
        phonopy.displacements, [[0, 0.0, np.sqrt(2) / 2, np.sqrt(2) / 2]]
    )

    phonopy = get_phonopy(atoms, symprec=1e-8)
    assert phonopy.symmetry.tolerance == 1e-8


def test_get_supercell_by_phonopy():
    atoms = bulk("Cu")
    get_atoms_supercell_by_phonopy(atoms, np.eye(3) * 2)

    cell = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]

    supercell = get_atoms_supercell_by_phonopy(atoms, cell)
    assert_almost_equal(np.diag(np.diag(supercell.cell)), supercell.cell)
