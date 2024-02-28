import pytest

ct = pytest.importorskip("covalent")

from ase.build import bulk

# from quacc import flow
# from quacc.recipes.emt.core import relax_job
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


def test_functools(tmp_path, monkeypatch):
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


# def test_copy_files(tmp_path, monkeypatch):
#     monkeypatch.chdir(tmp_path)
#     atoms = bulk("Cu")

#     @flow
#     def myflow(atoms):
#         result1 = relax_job(atoms)
#         return relax_job(result1["atoms"], copy_files={result1["dir_name"]: "opt.*"})

#     dispatch_id = ct.dispatch(myflow)(atoms)
#     output = ct.get_result(dispatch_id, wait=True)
#     assert "atoms" in output.result


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
