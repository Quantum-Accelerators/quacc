from __future__ import annotations

import pytest

redun = pytest.importorskip("redun")
from ase.build import bulk


@pytest.fixture()
def scheduler():
    return redun.Scheduler()


from quacc import flow, job
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


@pytest.mark.parametrize("job_decorators", [None, {"relax_job": job()}])
def test_functools(tmp_path, monkeypatch, scheduler, job_decorators):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    result = scheduler.run(
        bulk_to_slabs_flow(
            atoms,
            run_static=False,
            job_params={"relax_job": {"opt_params": {"fmax": 0.1}}},
            job_decorators=job_decorators,
        )
    )
    assert len(result) == 4
    assert "atoms" in result[-1]
    assert result[-1]["parameters_opt"]["fmax"] == 0.1


def test_copy_files(tmp_path, monkeypatch, scheduler):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")

    @flow
    def myflow(atoms):
        result1 = relax_job(atoms)
        return relax_job(result1["atoms"], copy_files={result1["dir_name"]: "opt.*"})

    assert "atoms" in scheduler.run(myflow(atoms))


def test_phonon_flow(tmp_path, monkeypatch, scheduler):
    pytest.importorskip("phonopy")
    pytest.importorskip("seekpath")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = scheduler.run(phonon_flow(atoms))
    assert output["results"]["thermal_properties"]["temperatures"].shape == (101,)
