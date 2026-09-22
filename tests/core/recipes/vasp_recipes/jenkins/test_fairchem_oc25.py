"""Opt-in OC25 VASP smoke test using locally licensed potentials."""

from __future__ import annotations

import os

import numpy as np
import pytest
from ase.build import add_adsorbate, fcc111, molecule
from ase.constraints import FixAtoms
from ase.io import read

from quacc.recipes.vasp.fairchem import has_fairchem_oc, oc25_static_job

pytestmark = [
    pytest.mark.skipif(
        os.environ.get("QUACC_RUN_OC25_VASP") != "1",
        reason="Set QUACC_RUN_OC25_VASP=1 with a configured VASP installation.",
    ),
    pytest.mark.skipif(not has_fairchem_oc, reason="fairchem-data-oc not installed"),
]


def test_fairchem_oc25_static():
    if structure := os.environ.get("QUACC_OC25_STRUCTURE"):
        atoms = read(structure)
    else:
        atoms = fcc111("Cu", size=(3, 3, 3), vacuum=8.0)
        bottom = atoms.get_tags() == 3
        atoms.set_constraint(FixAtoms(mask=bottom))
        atoms.set_tags(np.where(bottom, 0, 1))
        water = molecule("H2O")
        water.set_tags(2)
        add_adsorbate(atoms, water, height=2.5, position="ontop")
        atoms.pbc = True

    output = oc25_static_job(atoms)
    assert output["state"] == "successful"
    assert output["name"] == "OC25 Static"
    assert output["parameters"]["nsw"] == 0
    assert output["parameters"]["ediff"] == 1e-6
    assert output["parameters"]["pp_version"] == "64"
    assert np.isfinite(output["results"]["energy"])
    forces = np.asarray(output["results"]["forces"])
    assert forces.shape == (len(atoms), 3)
    assert np.isfinite(forces).all()
    np.testing.assert_allclose(output["atoms"].positions, atoms.positions, atol=1e-10)
    np.testing.assert_allclose(output["atoms"].cell, atoms.cell, atol=1e-10)
    np.testing.assert_array_equal(output["atoms"].get_tags(), atoms.get_tags())
