import os

import pytest
from custodian import Custodian

from quacc.custodian.vasp import run_custodian
from quacc.defaults import custodian_settings


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


def test_run_vasp_custodian(monkeypatch):
    monkeypatch.setenv(
        "VASP_CUSTODIAN_SETTINGS",
        os.path.join(os.path.dirname(custodian_settings.__file__),
                     "vasp_custodian_settings.yaml"),
    )
    monkeypatch.setenv("VASP_PARALLEL_CMD", "mpirun")
    run_custodian()


def test_failed_custodian(monkeypatch):
    with pytest.raises(OSError):
        run_custodian()

    monkeypatch.setenv(
        "VASP_CUSTODIAN_SETTINGS",
        os.path.join(os.path.dirname(custodian_settings.__file__),
                     "vasp_custodian_settings.yaml"),
    )
    with pytest.raises(OSError):
        run_custodian()
