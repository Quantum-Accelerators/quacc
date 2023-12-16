from functools import partial

import pytest
from ase.build import bulk

from quacc.recipes.tblite.core import static_job as static_job_
from quacc.recipes.tblite.phonons import phonon_flow

pytest.importorskip("tblite.ase")
pytest.importorskip("phonopy")


def test_phonon_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(atoms, static_job=partial(static_job_, method="GFN1-xTB"))
    assert output["results"]["force_constants"].shape == (8, 8, 3, 3)
    assert len(output["results"]["thermal_properties"]["temperatures"]) == 101
