import pytest
from ase.build import bulk

pytest.importorskip("phonopy")

from quacc.recipes.emt.phonons import phonon_flow


def test_bulk_to_defects_flow(tmpdir):
    tmpdir.chdir()
    atoms = bulk("Cu")
    output = phonon_flow(atoms)
    assert output["thermal_properties"]["temperatures"].shape == (101,)
