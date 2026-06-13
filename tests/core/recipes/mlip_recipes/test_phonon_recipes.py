from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("matcalc")
pytest.importorskip("matgl")
pytest.importorskip("phonopy")
pytest.importorskip("seekpath")
from ase.build import bulk

from quacc.recipes.mlip.phonons import phonon_flow


def test_phonon_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(
        atoms,
        min_lengths=5.0,
        job_params={
            "all": {"library": "matcalc", "name": "TensorNet-PES-MatPES-PBE-2025.2",
        },
    )
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert len(output["results"]["thermal_properties"]["temperatures"]) == 101
