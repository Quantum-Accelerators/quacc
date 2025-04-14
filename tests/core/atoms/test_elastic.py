from __future__ import annotations

from pathlib import Path

import pytest
from ase.io import read
from numpy.testing import assert_equal

from quacc.atoms.deformation import make_deformations_from_bulk

FILE_DIR = Path(__file__).parent


def test_make_deformations_from_bulk():
    atoms = read(FILE_DIR / "ZnTe.cif.gz")
    atoms.info["test"] = "hi"
    deformations = make_deformations_from_bulk(atoms)
    assert len(deformations) == 24
    assert deformations[0].to_ase_atoms().get_volume() != pytest.approx(
        atoms.get_volume()
    )
    for deformation in deformations:
        assert_equal(
            deformation.to_ase_atoms().get_atomic_numbers(),
            [30, 30, 30, 30, 52, 52, 52, 52],
        )
        assert_equal(deformation.to_ase_atoms().get_chemical_formula(), "Te4Zn4")
        assert deformation.to_ase_atoms().info["test"] == "hi"
