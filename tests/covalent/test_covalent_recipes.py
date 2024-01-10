import pytest

ct = pytest.importorskip("covalent")

from ase.build import bulk

from quacc.recipes.emt.core import relax_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


def test_covalent_functools(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(bulk_to_slabs_flow)(
        atoms, job_params={"relax_job": {"opt_params": {"fmax": 0.1}}}, run_static=False
    )
    output = ct.get_result(dispatch_id, wait=True)
    assert output.status == "COMPLETED"
    assert len(output.result) == 4
    assert "atoms" in output.result[-1]
    assert output.result[-1]["fmax"] == 0.1


def test_phonon_flow(tmp_path, monkeypatch):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(phonon_flow)(atoms)
    output = ct.get_result(dispatch_id, wait=True)
    assert output.status == "COMPLETED"
    assert output.result["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )
