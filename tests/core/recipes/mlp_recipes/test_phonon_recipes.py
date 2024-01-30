import pytest

pytest.importorskip("torch")
pytest.importorskip("mace")
pytest.importorskip("phonopy")
from ase.build import bulk

from quacc.recipes.mlp.phonons import phonon_flow


def test_phonon_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(atoms, method="mace", min_length=5.0)
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert len(output["results"]["thermal_properties"]["temperatures"]) == 101
