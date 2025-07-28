from __future__ import annotations

import pytest

pytest.importorskip("pymatgen.analysis.defects")
pytest.importorskip("shakenbreak")

from ase.build import bulk

from quacc.recipes.emt.defects import bulk_to_defects_flow


def test_bulk_to_defects_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = bulk_to_defects_flow(
        atoms, job_params={"relax_job": {"opt_params": {"fmax": 5}}}
    )
    assert len(output) == 2
    assert len(output[0]["atoms"]) == 107

    atoms = bulk("Cu")
    output = bulk_to_defects_flow(
        atoms, job_params={"relax_job": {"opt_params": {"fmax": 5}}}, run_static=False
    )
    assert len(output) == 2
    assert len(output[0]["atoms"]) == 107
