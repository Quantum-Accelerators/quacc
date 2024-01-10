import pytest

prefect = pytest.importorskip("prefect")

from ase.build import bulk
from prefect.testing.utilities import prefect_test_harness

from quacc.recipes.emt.core import relax_job  # skipcq: PYL-C0412
from quacc.recipes.emt.slabs import bulk_to_slabs_flow  # skipcq: PYL-C0412


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    with prefect_test_harness():
        yield


def test_prefect_functools(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    result = bulk_to_slabs_flow(
        atoms, job_params={"relax_job": {"opt_params": {"fmax": 0.1}}}, run_static=False
    ).result()
    assert len(result) == 4
    assert [r.is_completed() for r in result]


def test_phonon_flow(tmp_path, monkeypatch):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    output = phonon_flow(atoms)
    assert output.result().is_completed()


def test_phonon_flow_multistep(tmp_path, monkeypatch):
    pytest.importorskip("phonopy")
    from quacc.recipes.emt.phonons import phonon_flow

    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    relaxed = relax_job(atoms)
    output = phonon_flow(relaxed["atoms"])
    output.result().is_completed()
