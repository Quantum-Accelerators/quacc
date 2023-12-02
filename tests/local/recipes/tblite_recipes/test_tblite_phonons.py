import os

import pytest
from ase.build import bulk

from quacc.recipes.tblite.phonons import phonon_flow

pytest.importorskip("tblite.ase")
pytest.importorskip("phonopy")


def test_phonon_flow(tmp_path, monkeypatch:
    monkeypatch.chdir(tmp_path, monkeypatch
    atoms = bulk("Cu")
    output = phonon_flow(atoms, static_job_kwargs={"method": "GFN1-xTB"})
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert len(output["results"]["thermal_properties"]["temperatures"]) == 101
