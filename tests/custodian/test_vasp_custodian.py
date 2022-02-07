import os
import pytest

from custodian import Custodian

from quacc.custodian.vasp import run_custodian
from quacc.defaults import custodian_settings


class MockRun:
    @staticmethod
    def run():
        return True


def mock_custodian_run(*args, **kwargs):
    # Instead of running the VASP-specific summarize_run(), we mock it with the
    # general calculator schema which does not require VASP files to be
    # in the working directory and will work with pytest.

    return MockRun()


@pytest.fixture(autouse=True)
def patch_custodian_run(monkeypatch):
    # Monkeypatch the summarize_run() function so that we aren't relying on real
    # VASP files to be in the working directory during the test. Note that even though
    # summarize_run() is a function in the quacc.schemas.vasp module, we modify it
    # only in quacc.recipes.vasp.core/.slabs because otherwise it will not work properly.
    monkeypatch.setattr(Custodian, "run", mock_custodian_run)


def test_run_vasp_custodian(monkeypatch):

    monkeypatch.setenv(
        "VASP_CUSTODIAN_SETTINGS",
        os.path.join(
            os.path.dirname(custodian_settings.__file__), "vasp_custodian_settings.yaml"
        ),
    )
    monkeypatch.setenv("VASP_PARALLEL_CMD", "mpirun")
    run_custodian()


def test_failed_custodian(monkeypatch):
    with pytest.raises(OSError):
        run_custodian()

    monkeypatch.setenv(
        "VASP_CUSTODIAN_SETTINGS",
        os.path.join(
            os.path.dirname(custodian_settings.__file__), "vasp_custodian_settings.yaml"
        ),
    )
    with pytest.raises(OSError):
        run_custodian()
