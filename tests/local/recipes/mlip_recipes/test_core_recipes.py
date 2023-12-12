import numpy as np
import pytest
from ase.build import bulk

from quacc.recipes.mlip.core import relax_job, static_job


@pytest.mark.parametrize("method", ["chgnet", "m3gnet", "umace"])
@pytest.mark.parametrize(
    "energy", [-4.083308219909668, np.array([-4.0938973]), -4.083906650543213]
)
def test_relax_job(tmp_path, monkeypatch, method, energy):
    monkeypatch.chdir(tmp_path)
    if method == "chgnet":
        pytestmark = pytest.importorskip("chgnet")
    elif method == "m3gnet":
        pytestmark = pytest.importorskip("matgl")
        pytestmark = pytest.importorskip("matcalc")
    elif method == "umace":
        pytestmark = pytest.importorskip("mace")
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = static_job(atoms, method=method)
    assert output["results"]["energy"] == pytest.approx(energy)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert np.shape(output["results"]["stress"]) == (3, 3)
    assert output["atoms"] == atoms


@pytest.mark.parametrize("method", ["chgnet", "m3gnet", "umace"])
@pytest.mark.parametrize(
    "energy", [-32.662731, np.array([-32.747219]), -32.66752624511719]
)
def test_relax_job(tmp_path, monkeypatch, method, energy):
    monkeypatch.chdir(tmp_path)
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
    assert output["results"]["energy"] == pytest.approx(energy)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert np.shape(output["results"]["stress"]) == (3, 3)
    assert output["atoms"] != atoms
