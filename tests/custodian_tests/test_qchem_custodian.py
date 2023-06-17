import pytest
from custodian import Custodian

try:
    import openbabel as ob
except ImportError:
    ob = None


class MockRun:
    # Add a mock Custodian.run() function

    @staticmethod
    def run():
        return True


def mock_custodian_run(*args, **kwargs):
    # Instead of running Custodian, we will mock it to return True
    # when .run() is called
    return MockRun()


@pytest.fixture(autouse=True)
def patch_custodian_run(monkeypatch):
    # Monkeypatch the Custodian.run() function so that it doesn't actually
    # launch Custodian during a test
    monkeypatch.setattr(Custodian, "run", mock_custodian_run)


@pytest.mark.skipif(
    ob is None,
    reason="Openbabel needed for test.",
)
def test_run_qchem_custodian(monkeypatch):
    from quacc.custodian.qchem import run_custodian

    run_custodian()

    run_custodian(
        qchem_cores=40,
        qchem_cmd="qchem -save",
        qchem_local_scratch="/not_tmp",
        qchem_custodian_max_errors=20,
    )

    run_custodian(qchem_use_error_handlers=False)
