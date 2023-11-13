import pytest
from custodian import Custodian

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

    with pytest.raises(ValueError):
        run_custodian(vasp_custodian_handlers="cow")

    with pytest.raises(ValueError):
        run_custodian(vasp_custodian_validators=["cow"])
