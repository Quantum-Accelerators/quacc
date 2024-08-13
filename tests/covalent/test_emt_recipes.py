from __future__ import annotations

import pytest

ct = pytest.importorskip("covalent")

from ase.build import bulk

from quacc import flow, job
from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


@pytest.mark.parametrize("job_decorators", [None, {"relax_job": job()}])
def test_functools(tmp_path, monkeypatch, job_decorators):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(bulk_to_slabs_flow)(
        atoms,
        run_static=False,
        job_params={"relax_job": {"opt_params": {"fmax": 0.1}}},
        job_decorators=job_decorators,
    )
    output = ct.get_result(dispatch_id, wait=True)
    assert output.status == "COMPLETED"
    assert len(output.result) == 4
    assert "atoms" in output.result[-1]
    assert output.result[-1]["parameters_opt"]["fmax"] == 0.1


def test_copy_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")

    @flow
    def myflow(atoms):
        result1 = relax_job(atoms)
        return relax_job(result1["atoms"], copy_files={result1["dir_name"]: "opt.*"})

    dispatch_id = ct.dispatch(myflow)(atoms)
    output = ct.get_result(dispatch_id, wait=True)
    assert "atoms" in output.result


def test_phonon_flow(tmp_path, monkeypatch):
    pytest.importorskip("phonopy")
    pytest.importorskip("seekpath")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(phonon_flow)(atoms)
    output = ct.get_result(dispatch_id, wait=True)
    assert output.status == "COMPLETED"
    assert output.result["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )
