import pytest
from ase.build import bulk

from quacc import SETTINGS
from quacc.recipes.emt.core import relax_job

ct = pytest.importorskip("covalent")
pytestmark = pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE != "covalent",
    reason="This test requires the Covalent workflow engine",
)


def test_phonon_flow(tmp_path):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    tmp_path.chdir()
    atoms = bulk("Cu")
    dispatch_id = ct.dispatch(phonon_flow)(atoms)
    output = ct.get_result(dispatch_id, wait=True)
    assert output.status == "COMPLETED"
    assert output.result["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )


def test_phonon_flow_multistep(tmp_path):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    tmp_path.chdir()
    atoms = bulk("Cu")
    relaxed = relax_job(atoms)
    dispatch_id = ct.dispatch(phonon_flow)(relaxed["atoms"])
    output = ct.get_result(dispatch_id, wait=True)
    assert output.status == "COMPLETED"
    assert output.result["results"]["thermal_properties"]["temperatures"].shape == (
        101,
    )
