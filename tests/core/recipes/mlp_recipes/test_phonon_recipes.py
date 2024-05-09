from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("matgl")
pytest.importorskip("phonopy")
import numpy as np
from ase.build import bulk

from quacc.recipes.mlp.phonons import phonon_flow


def _set_dtype(size, type_="float"):
    globals()[f"{type_}_th"] = getattr(torch, f"{type_}{size}")
    globals()[f"{type_}_np"] = getattr(np, f"{type_}{size}")
    torch.set_default_dtype(getattr(torch, f"float{size}"))


@pytest.mark.xfail(strict=False)
def test_phonon_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _set_dtype(64)
    atoms = bulk("Cu")
    output = phonon_flow(atoms, method="m3gnet", min_lengths=5.0)
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert len(output["results"]["thermal_properties"]["temperatures"]) == 101


@pytest.mark.xfail(strict=False)
def test_phonon_flow_dispersion(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _set_dtype(64)
    atoms = bulk("Cu")
    output = phonon_flow(
        atoms,
        method="m3gnet",
        min_lengths=5.0,
        job_params={"all": {"compute_stress": False}},
    )
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert len(output["results"]["thermal_properties"]["temperatures"]) == 101
