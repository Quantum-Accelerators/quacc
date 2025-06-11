from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")

from importlib.util import find_spec

import numpy as np
from ase.build import bulk

from quacc.recipes.mlp.core import relax_job, static_job

methods = []
if has_mace := find_spec("mace"):
    methods.append("mace-mp")

if has_matgl := find_spec("matgl"):
    methods.extend(["m3gnet", "chgnet", "tensornet"])

if has_sevennet := find_spec("sevenn"):
    methods.append("sevennet")

if has_orb := find_spec("orb_models"):
    methods.append("orb")

if find_spec("fairchem"):
    from huggingface_hub.utils._auth import get_token

    if get_token():
        methods.append("fairchem")


def _set_dtype(size, type_="float"):
    globals()[f"{type_}_th"] = getattr(torch, f"{type_}{size}")
    globals()[f"{type_}_np"] = getattr(np, f"{type_}{size}")
    torch.set_default_dtype(getattr(torch, f"float{size}"))


@pytest.mark.parametrize("method", methods)
def test_static_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp":
        _set_dtype(64)
    else:
        _set_dtype(32)

    if method == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}
    else:
        calc_kwargs = {}

    ref_energy = {
        "chgnet": -3.7441039085388184,
        "m3gnet": -3.7398147583007812,
        "tensornet": 5.0,
        "mace-mp": -4.097862720291976,
        "sevennet": -4.096191883087158,
        "orb": -4.093477725982666,
        "fairchem": -3.7579006783217954,
    }
    atoms = bulk("Cu")
    output = static_job(atoms, method=method, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (1, 3)
    assert output["atoms"] == atoms


@pytest.mark.skipif(has_sevennet is None, reason="sevennet not installed")
def test_static_job_with_dict_kwargs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    # Make sure that pick_calculator works even with dictionary kwargs
    static_job(atoms, method="sevennet", sevennet_config={"test": 1})


def test_relax_job_missing_pynanoflann(monkeypatch):
    def mock_find_spec(name):
        if name == "pynanoflann":
            return None
        return find_spec(name)

    import quacc.recipes.mlp._base

    quacc.recipes.mlp._base.pick_calculator.__wrapped__.cache_clear()
    monkeypatch.setattr("importlib.util.find_spec", mock_find_spec)
    monkeypatch.setattr("quacc.recipes.mlp._base.find_spec", mock_find_spec)
    with pytest.raises(ImportError, match=r"orb-models requires pynanoflann"):
        relax_job(bulk("Cu"), method="orb")


@pytest.mark.parametrize("method", methods)
def test_relax_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp":
        _set_dtype(64)
    else:
        _set_dtype(32)

    if method == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}
    else:
        calc_kwargs = {}

    ref_energy = {
        "chgnet": -29.952457427978516,
        "m3gnet": -29.9184513092041,
        "mace-mp": -32.78264569638644,
        "tensornet": 5.0,
        "sevennet": -32.76924133300781,
        "orb": -32.7361946105957,
        "fairchem": -30.004380887389797,
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
    output = relax_job(atoms, method="mace-mp", dispersion=True)
    assert output["results"]["energy"] == pytest.approx(-37.4518034464096)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())


@pytest.mark.parametrize("method", methods)
def test_relax_cell_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace-mp":
        _set_dtype(64)
    else:
        _set_dtype(32)

    if method == "fairchem":
        # Note that for this to work, you need HF_TOKEN env variable set!
        calc_kwargs = {"name_or_path": "uma-s-1", "task_name": "omat"}
    else:
        calc_kwargs = {}

    ref_energy = {
        "chgnet": -29.966711044311523,
        "m3gnet": -29.933645248413086,
        "mace-mp": -32.8069374165035,
        "tensornet": 5.0,
        "sevennet": -32.76963806152344,
        "orb": -32.73428726196289,
        "fairchem": -30.005004590392726,
    }

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method, relax_cell=True, **calc_kwargs)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method], rel=1e-4)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() != pytest.approx(atoms.get_volume())
