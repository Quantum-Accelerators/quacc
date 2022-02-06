import pytest


@pytest.fixture(autouse=True)
def set_env(monkeypatch):
    # Set environment variables in pytest so the tests can run
    monkeypatch.setenv("VASP_PARALLEL_CMD", "")
    monkeypatch.setenv("VASP_PP_PATH", ".")
    