from __future__ import annotations

import pytest
from custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler

from quacc.calculators.vasp import vasp_custodian
from quacc.calculators.vasp.vasp_custodian import run_custodian


def mock_custodian_run(*args, **kwargs):
    # Instead of running Custodian, we will mock it to return True
    # when .run() is called

    class MockRun:
        # Add a mock Custodian.run() function

        @staticmethod
        def run():
            return True

    return MockRun()


@pytest.fixture(autouse=True)
def patch_custodian_run(monkeypatch):
    # Monkeypatch the Custodian.run() function so that it doesn't actually
    # launch Custodian during a test

    monkeypatch.setattr(Custodian, "run", mock_custodian_run)


def test_run_vasp_custodian(monkeypatch):
    monkeypatch.setenv("VASP_PARALLEL_CMD", "fake-mpirun")
    run_custodian()

    run_custodian(vasp_custodian_wall_time=1)
    run_custodian(vasp_custodian_handlers=None, vasp_custodian_validators=None)

    with pytest.raises(ValueError, match="Unknown VASP error handler"):
        run_custodian(vasp_custodian_handlers="cow")

    with pytest.raises(ValueError, match="Unknown VASP validator"):
        run_custodian(vasp_custodian_validators=["cow"])


def test_vasp_error_exclude(monkeypatch):
    original_error_msgs = VaspErrorHandler.error_msgs
    seen_handlers = []

    class MockVaspErrorHandler:
        error_msgs = original_error_msgs
        is_monitor = False

        def __init__(self, **kwargs):
            self.kwargs = kwargs
            seen_handlers.append(self)

    monkeypatch.setattr(vasp_custodian, "VaspErrorHandler", MockVaspErrorHandler)

    run_custodian(vasp_custodian_vasp_error_exclude=["bravais"])

    handler = seen_handlers[0]
    assert "bravais" not in handler.kwargs["errors_subset_to_catch"]
    assert set(handler.kwargs["errors_subset_to_catch"]) == set(original_error_msgs) - {
        "bravais"
    }

    with pytest.raises(ValueError, match="Unknown VASP errors to exclude"):
        run_custodian(vasp_custodian_vasp_error_exclude=["not-a-vasp-error"])
