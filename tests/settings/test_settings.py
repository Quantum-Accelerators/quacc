import os

import pytest

from quacc import SETTINGS


@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions",
)
def test_settings():
    # These settings were defined in the test.yaml GitHub actions workflow

    assert SETTINGS.SCRATCH_DIR == "/tmp"
    assert SETTINGS.VASP_VERBOSE is False
    assert SETTINGS.VASP_CUSTODIAN is True
