import numpy as np
import pytest
from ase.build import bulk

from quacc.recipes.mlip.core import relax_job, static_job


@pytest.mark.parametrize("method", ["chgnet", "m3gnet", "umace"])
def test_static_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)
    ref_energy = {
        "chgnet": -4.083308219909668,
        "m3gnet": -4.0938973,
        "umace": -4.083906650543213,
    }
    if method == "chgnet":
        pytestmark = pytest.importorskip("chgnet")
    elif method == "m3gnet":
        pytestmark = pytest.importorskip("matgl")
        pytestmark = pytest.importorskip("matcalc")
    elif method == "umace":
        pytestmark = pytest.importorskip("mace")
    atoms = bulk("Cu")
    output = static_job(atoms, method=method)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (1, 3)
    assert output["atoms"] == atoms


@pytest.mark.parametrize("method", ["chgnet", "m3gnet", "umace"])
def test_relax_job(tmp_path, monkeypatch, method):
    monkeypatch.chdir(tmp_path)
    ref_energy = {
        "chgnet": -32.662731,
        "m3gnet": -32.747219,
        "umace": -32.66752624511719,
    }
    if method == "chgnet":
        pytestmark = pytest.importorskip("chgnet")
    elif method == "m3gnet":
        pytestmark = pytest.importorskip("matgl")
        pytestmark = pytest.importorskip("matcalc")
    elif method == "umace":
        pytestmark = pytest.importorskip("mace")

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(atoms, method=method)
    assert output["results"]["energy"] == pytest.approx(ref_energy[method])
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
