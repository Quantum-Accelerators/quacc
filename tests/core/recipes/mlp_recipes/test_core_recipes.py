import numpy as np
import pytest
from ase.build import bulk

from quacc.recipes.mlp.core import relax_job, static_job

torch = pytest.importorskip("torch")


def test_bad_method():
    atoms = bulk("Cu")
    pytestmark = pytest.importorskip("chgnet")
    with pytest.raises(ValueError):
        static_job(atoms, method="bad_method")


def set_default_dtype(method):
    if method != "mace":
        torch.set_default_dtype(getattr(torch, "float32", "float32"))


@pytest.mark.parametrize("method", ["chgnet", "m3gnet", "mace"])
def test_static_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)
    set_default_dtype(method)

    ref_energy = {
        "chgnet": -4.083308219909668,
        "m3gnet": -4.0938973,
        "mace": -4.083906650543213,
    }
    if method == "chgnet":
        pytestmark = pytest.importorskip("chgnet")
    elif method == "m3gnet":
        pytestmark = pytest.importorskip("matgl")
    elif method == "mace":
        pytestmark = pytest.importorskip("mace")
    atoms = bulk("Cu")
    output = static_job(atoms, method=method)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (1, 3)
    assert output["atoms"] == atoms


@pytest.mark.parametrize("method", ["chgnet", "m3gnet", "mace"])
def test_relax_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)
    set_default_dtype(method)

    ref_energy = {
        "chgnet": -32.665626525878906,
        "m3gnet": -32.749088287353516,
        "mace": -32.670471191406259,
    }
    if method == "chgnet":
        pytestmark = pytest.importorskip("chgnet")
    elif method == "m3gnet":
        pytestmark = pytest.importorskip("matgl")
    elif method == "mace":
        pytestmark = pytest.importorskip("mace")

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())


@pytest.mark.parametrize("method", ["chgnet", "m3gnet", "mace"])
def test_relax_cell_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)
    set_default_dtype(method)

    ref_energy = {
        "chgnet": -32.6676139831543,
        "m3gnet": -32.74995040893555,
        "mace": -32.67771911621094,
    }
    if method == "chgnet":
        pytestmark = pytest.importorskip("chgnet")
    elif method == "m3gnet":
        pytestmark = pytest.importorskip("matgl")
    elif method == "mace":
        pytestmark = pytest.importorskip("mace")

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method, relax_cell=True)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() != pytest.approx(atoms.get_volume())
