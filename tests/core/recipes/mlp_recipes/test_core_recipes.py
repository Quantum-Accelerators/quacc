import pytest

torch = pytest.importorskip("torch")

import numpy as np
from ase.build import bulk

from quacc.recipes.mlp.core import relax_job, static_job

methods = []
try:
    import mace

    methods.append("mace")

except ImportError:
    mace = None
# try:
#     import matgl

#     methods.append("m3gnet")
# except ImportError:
#     matgl = None
try:
    import chgnet

    methods.append("chgnet")
except ImportError:
    chgnet = None


@pytest.mark.skipif(chgnet is None, reason="chgnet not installed")
def test_bad_method():
    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        static_job(atoms, method="bad_method")


def _set_dtype(size, type_="float"):
    globals()[f"{type_}_th"] = getattr(torch, f"{type_}{size}")
    globals()[f"{type_}_np"] = getattr(np, f"{type_}{size}")
    torch.set_default_dtype(getattr(torch, f"float{size}"))


@pytest.mark.parametrize("method", methods)
def test_static_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace":
        _set_dtype(64)
    else:
        _set_dtype(32)

    ref_energy = {
        "chgnet": -4.083308219909668,
        "m3gnet": -4.0938973,
        "mace": -4.083906650543213,
    }
    atoms = bulk("Cu")
    output = static_job(atoms, method=method)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (1, 3)
    assert output["atoms"] == atoms


@pytest.mark.parametrize("method", methods)
def test_relax_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace":
        _set_dtype(64)
    else:
        _set_dtype(32)
    ref_energy = {
        "chgnet": -32.665626525878906,
        "m3gnet": -32.749088287353516,
        "mace": -32.670471191406259,
    }

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())


def test_relax_job_dispersion(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    _set_dtype(64)

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method="mace", dispersion=True)
    assert output["results"]["energy"] == pytest.approx(-37.33948477096204)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())


@pytest.mark.parametrize("method", methods)
def test_relax_cell_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)

    if method == "mace":
        _set_dtype(64)
    else:
        _set_dtype(32)

    ref_energy = {
        "chgnet": -32.6676139831543,
        "m3gnet": -32.74995040893555,
        "mace": -32.67771911621094,
    }

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method, relax_cell=True)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() != pytest.approx(atoms.get_volume())
