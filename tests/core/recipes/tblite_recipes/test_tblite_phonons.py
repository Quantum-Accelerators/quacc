import pytest

pytest.importorskip("tblite.ase")
pytest.importorskip("phonopy")
from ase.build import bulk

from quacc.recipes.tblite.phonons import phonon_flow


def test_phonon_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(
        atoms, min_lengths=5.0, job_params={"static_job": {"method": "GFN1-xTB"}}
    )
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert len(output["results"]["thermal_properties"]["temperatures"]) == 101
