import pytest

pytest.importorskip("phonopy")
from ase.build import bulk

from quacc.recipes.emt.phonons import phonon_flow


def test_phonon_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(atoms, min_length=5.0)
    assert output["results"]["thermal_properties"]["temperatures"].shape == (101,)
    assert output["results"]["thermal_properties"]["temperatures"][0] == 0
    assert output["results"]["thermal_properties"]["temperatures"][-1] == 1000
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert "mesh_properties" in output["results"]

    output = phonon_flow(atoms, min_length=None, t_min=10, t_max=20, t_step=5)
    assert output["results"]["thermal_properties"]["temperatures"].shape == (3,)
    assert output["results"]["thermal_properties"]["temperatures"][0] == 10
    assert output["results"]["thermal_properties"]["temperatures"][-1] == 20
    assert output["results"]["force_constants"].shape == (1, 1, 3, 3)
    assert "mesh_properties" in output["results"]
