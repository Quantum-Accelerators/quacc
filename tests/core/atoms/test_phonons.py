from __future__ import annotations

import pytest

pytest.importorskip("phonopy")
pytest.importorskip("seekpath")

import numpy as np
from ase.build import bulk
from ase.constraints import FixAtoms
from numpy.testing import assert_almost_equal, assert_array_equal

from quacc.atoms.phonons import get_phonopy


def test_get_phonopy():
    atoms = bulk("Cu")
    phonopy, _ = get_phonopy(atoms)
    assert_array_equal(phonopy.supercell_matrix, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    phonopy, _ = get_phonopy(atoms, min_lengths=5)
    assert_array_equal(phonopy.supercell_matrix, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])

    phonopy, _ = get_phonopy(atoms, min_lengths=[5, 10, 5])
    assert_array_equal(phonopy.supercell_matrix, [[2, 0, 0], [0, 4, 0], [0, 0, 2]])

    phonopy, _ = get_phonopy(atoms, displacement=1)
    assert_almost_equal(
        phonopy.displacements, [[0, 0.0, np.sqrt(2) / 2, np.sqrt(2) / 2]]
    )

    phonopy, _ = get_phonopy(atoms, symprec=1e-8)
    assert phonopy.symmetry.tolerance == 1e-8

    atoms = bulk("Cu") * (2, 2, 2)

    atoms.set_constraint(FixAtoms(indices=[0, 1, 2, 3]))

    phonopy, fixed_atoms = get_phonopy(atoms, min_lengths=5)

    assert len(fixed_atoms) == 4
