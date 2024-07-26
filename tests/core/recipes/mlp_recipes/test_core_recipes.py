from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from importlib.util import find_spec

import numpy as np
from ase.build import bulk

from quacc.recipes.mlp.core import relax_job, static_job

methods = []
if has_mace := find_spec("mace"):
    methods.append("mace-mp-0")

if has_matgl := find_spec("matgl"):
    methods.append("m3gnet")

if has_chgnet := find_spec("chgnet"):
    methods.append("chgnet")


@pytest.mark.skipif(has_chgnet is None, reason="chgnet not installed")
def test_bad_method():
    atoms = bulk("Cu")
    with pytest.raises(ValueError, match="Unrecognized method='bad_method'"):
        static_job(atoms, method="bad_method")


def _set_dtype(size, type_="float"):
    globals()[f"{type_}_th"] = getattr(torch, f"{type_}{size}")
    globals()[f"{type_}_np"] = getattr(np, f"{type_}{size}")
    torch.set_default_dtype(getattr(torch, f"float{size}"))


@pytest.mark.parametrize("method", methods)
def test_static_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp-0":
        _set_dtype(64)
    else:
        _set_dtype(32)

    ref_energy = {
        "chgnet": -4.083308219909668,
        "m3gnet": -4.0938973,
        "mace-mp-0": -4.083906650543213,
    }
    atoms = bulk("Cu")
    output = static_job(atoms, method=method)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (1, 3)
    assert output["atoms"] == atoms


@pytest.mark.parametrize("method", methods)
def test_relax_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp-0":
        _set_dtype(64)
    else:
        _set_dtype(32)
    ref_energy = {
        "chgnet": -32.665428161621094,
        "m3gnet": -32.75003433227539,
        "mace-mp-0": -32.6711566550002,
    }

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())


@pytest.mark.skipif(has_mace is None, reason="Needs MACE")
@pytest.mark.importorskip("torch_dftd", reason="Needs torch-dftd")
def test_relax_job_dispersion(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    _set_dtype(64)

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method="mace-mp-0", dispersion=True)
    assert output["results"]["energy"] == pytest.approx(-37.340311589504076)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())


@pytest.mark.parametrize("method", methods)
def test_relax_cell_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp-0":
        _set_dtype(64)
    else:
        _set_dtype(32)

    ref_energy = {
        "chgnet": -32.66698455810547,
        "m3gnet": -32.750858306884766,
        "mace-mp-0": -32.67840391814377,
    }

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method, relax_cell=True)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() != pytest.approx(atoms.get_volume())
