"""Tests for optional file I/O in MLIP recipes."""

from __future__ import annotations

from ase.build import bulk
from ase.calculators.emt import EMT

from quacc.recipes.mlip.core import relax_job


def test_relax_job_without_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        "quacc.recipes.mlip.core.pick_calculator", lambda *_args, **_kwargs: EMT()
    )

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    output = relax_job(atoms, library="matcalc", write_files=False)

    assert len(output["trajectory"]) == 1
    assert output["trajectory"][0] == output["atoms"]
    assert not list(tmp_path.rglob("opt.log*"))
    assert not list(tmp_path.rglob("opt.json*"))
    assert not list(tmp_path.rglob("opt.traj*"))
