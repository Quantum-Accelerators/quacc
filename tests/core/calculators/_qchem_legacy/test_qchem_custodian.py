pytest.importorskip("openbabel")
import pytest
from custodian import Custodian

from quacc.calculators._qchem_legacy.qchem_custodian import run_custodian


def mock_custodian_run(*args, **kwargs):
    """Instead of running Custodian, we will mock it to return True when .run() is called."""

    class MockRun:
        "Mock Custodian run() function"

        @staticmethod
        def run():
            return True

    return MockRun()


@pytest.fixture(autouse=True)
def patch_custodian_run(monkeypatch):
    """Monkeypatch the Custodian.run() function so that it doesn't actually launch Custodian during a test."""

    monkeypatch.setattr(Custodian, "run", mock_custodian_run)


def test_run_qchem_custodian(monkeypatch):
    run_custodian()

    run_custodian(
        qchem_cores=40,
        qchem_cmd="qchem -save",
        qchem_local_scratch="/not_tmp",
        qchem_custodian_max_errors=20,
    )

    run_custodian(qchem_use_error_handlers=False)
