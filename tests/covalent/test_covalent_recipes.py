import pytest
from ase.build import bulk

ct = pytest.importorskip("covalent")


def test_phonon_flow(tmpdir):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    tmpdir.chdir()
    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(phonon_flow)(atoms)
    output = ct.get_result(dispatch_id, wait=True)
    assert output.status == "COMPLETED"
    assert output.result["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )
