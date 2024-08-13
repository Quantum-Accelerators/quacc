from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("phonopy")
pytest.importorskip("seekpath")

from ase.build import bulk

from quacc.recipes.emt.phonons import phonon_flow


def test_phonon_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(atoms, min_lengths=5.0)
    assert output["results"]["thermal_properties"]["temperatures"].shape == (101,)
    assert output["results"]["thermal_properties"]["temperatures"][0] == 0
    assert output["results"]["thermal_properties"]["temperatures"][-1] == 1000
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert "mesh_properties" in output["results"]
    assert "total_dos" in output["results"]
    results_dir = Path(output["dir_name"])
    assert Path(results_dir / "phonopy.yaml.gz").is_file()
    assert Path(results_dir, "phonopy_auto_band_structure.yaml.gz").is_file()

    atoms = bulk("Cu")
    output = phonon_flow(atoms, supercell_matrix=((2, 0, 0), (0, 2, 0), (0, 0, 2)))
    assert output["results"]["thermal_properties"]["temperatures"].shape == (101,)
    assert output["results"]["thermal_properties"]["temperatures"][0] == 0
    assert output["results"]["thermal_properties"]["temperatures"][-1] == 1000
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert "mesh_properties" in output["results"]
    assert "total_dos" in output["results"]


def test_phonon_flow_v2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu") * (2, 2, 2)
    output = phonon_flow(atoms, min_lengths=None, t_min=10, t_max=20, t_step=5)
    assert output["results"]["thermal_properties"]["temperatures"].shape == (3,)
    assert output["results"]["thermal_properties"]["temperatures"][0] == 10
    assert output["results"]["thermal_properties"]["temperatures"][-1] == 20
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert "mesh_properties" in output["results"]


def test_phonon_flow_v3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.2
    output = phonon_flow(atoms, min_lengths=5.0)
    assert output["results"]["thermal_properties"]["temperatures"].shape == (101,)
    assert output["results"]["thermal_properties"]["temperatures"][0] == 0
    assert output["results"]["thermal_properties"]["temperatures"][-1] == 1000
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert "mesh_properties" in output["results"]
    assert output["atoms"] == atoms
