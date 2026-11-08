from __future__ import annotations

from pathlib import Path

import pytest
from ase.io import read
from numpy.testing import assert_allclose, assert_equal

from quacc.atoms.deformation import (
    make_deformations_from_bulk,
    make_volume_strains_from_bulk,
)

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


def test_make_volume_strains_from_bulk():
    atoms = read(FILE_DIR / "ZnTe.cif.gz")
    atoms.info["test"] = "hi"
    volume_strains = [-0.075, -0.01, 0.01, 0.075]
    strained_atoms = make_volume_strains_from_bulk(atoms, volume_strains=volume_strains)
    assert len(strained_atoms) == len(volume_strains)
    for strain, strained in zip(volume_strains, strained_atoms, strict=True):
        assert strained.get_volume() == pytest.approx((1 + strain) * atoms.get_volume())
        # Isotropic: every lattice vector is scaled by the same factor
        assert_allclose(strained.get_cell(), atoms.get_cell() * (1 + strain) ** (1 / 3))
        assert_allclose(strained.get_scaled_positions(), atoms.get_scaled_positions())
        assert_equal(strained.get_chemical_formula(), "Te4Zn4")
        assert strained.info["test"] == "hi"

    assert len(make_volume_strains_from_bulk(atoms)) == 10

    with pytest.raises(ValueError, match="greater than -1"):
        make_volume_strains_from_bulk(atoms, volume_strains=[-1.0])
