from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("matgl")
pytest.importorskip("phonopy")
pytest.importorskip("seekpath")
from ase.build import bulk

from quacc.recipes.mlp.phonons import phonon_flow


def test_phonon_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(atoms, method="m3gnet", min_lengths=5.0)
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert len(output["results"]["thermal_properties"]["temperatures"]) == 101


def test_phonon_flow_dispersion(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(
        atoms,
        method="m3gnet",
        min_lengths=5.0,
        job_params={"all": {"compute_stress": False}},
    )
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert len(output["results"]["thermal_properties"]["temperatures"]) == 101
