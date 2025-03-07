from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from importlib.util import find_spec
from pathlib import Path

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

if has_sevennet := find_spec("sevenn"):
    methods.append("sevennet")

if has_orb := find_spec("orb_models"):
    methods.append("orb")

if has_fairchem := find_spec("fairchem"):
    methods.append("fairchem")


def _set_dtype(size, type_="float"):
    globals()[f"{type_}_th"] = getattr(torch, f"{type_}{size}")
    globals()[f"{type_}_np"] = getattr(np, f"{type_}{size}")
    torch.set_default_dtype(getattr(torch, f"float{size}"))


@pytest.mark.skipif(has_chgnet is None, reason="chgnet not installed")
def test_bad_method():
    atoms = bulk("Cu")
    with pytest.raises(ValueError, match="Unrecognized method='bad_method'"):
        static_job(atoms, method="bad_method")


@pytest.mark.parametrize("method", methods)
def test_static_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp-0":
        _set_dtype(64)
    else:
        _set_dtype(32)

    if method == "fairchem":
        calc_kwargs = {
            "checkpoint_path": Path(__file__).parent / "eqV2_31M_omat_mp_salex.pt"
        }
    else:
        calc_kwargs = {}

    ref_energy = {
        "chgnet": -4.083308219909668,
        "m3gnet": -4.0938973,
        "mace-mp-0": -4.097862720291976,
        "sevennet": -4.096191883087158,
        "orb": -4.093477725982666,
        "fairchem": -4.098316669464111,
    }
    atoms = bulk("Cu")
    output = static_job(atoms, method=method, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (1, 3)
    assert output["atoms"] == atoms


def test_relax_job_missing_pynanoflann(monkeypatch):
    def mock_find_spec(name):
        if name == "pynanoflann":
            return None
        return find_spec(name)

    import quacc.recipes.mlp._base

    quacc.recipes.mlp._base.pick_calculator.cache_clear()
    monkeypatch.setattr("importlib.util.find_spec", mock_find_spec)
    monkeypatch.setattr("quacc.recipes.mlp._base.find_spec", mock_find_spec)
    with pytest.raises(ImportError, match=r"orb-models requires pynanoflann"):
        relax_job(bulk("Cu"), method="orb")


@pytest.mark.parametrize("method", methods)
def test_relax_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp-0":
        _set_dtype(64)
    else:
        _set_dtype(32)

    if method == "fairchem":
        calc_kwargs = {
            "checkpoint_path": Path(__file__).parent / "eqV2_31M_omat_mp_salex.pt"
        }
    else:
        calc_kwargs = {}

    ref_energy = {
        "chgnet": -32.665428161621094,
        "m3gnet": -32.75003433227539,
        "mace-mp-0": -32.78264569638644,
        "sevennet": -32.76924133300781,
        "orb": -32.7361946105957,
        "fairchem": -32.80327224731445,
    }

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())


@pytest.mark.skipif(has_mace is None, reason="Needs MACE")
@pytest.mark.skipif(find_spec("torch_dftd") is None, reason="Needs torch-dftd")
def test_relax_job_dispersion(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    _set_dtype(64)

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method="mace-mp-0", dispersion=True)
    assert output["results"]["energy"] == pytest.approx(-37.4518034464096)
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

    if method == "fairchem":
        calc_kwargs = {
            "checkpoint_path": Path(__file__).parent / "eqV2_31M_omat_mp_salex.pt"
        }
    else:
        calc_kwargs = {}

    ref_energy = {
        "chgnet": -32.66698455810547,
        "m3gnet": -32.750858306884766,
        "mace-mp-0": -32.8069374165035,
        "sevennet": -32.76963806152344,
        "orb": -32.73428726196289,
        "fairchem": -32.82823944091797,
    }

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method, relax_cell=True, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() != pytest.approx(atoms.get_volume())
