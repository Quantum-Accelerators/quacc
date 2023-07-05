import pytest

from quacc import SETTINGS


@pytest.fixture(autouse=True)
def patch_settings(monkeypatch):
    # Monkeypatch the the check convergence setting
    monkeypatch.setattr(SETTINGS, "CHECK_CONVERGENCE", False)
