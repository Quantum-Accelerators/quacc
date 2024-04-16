from __future__ import annotations

import shutil
from pathlib import Path

import pytest

FILE_DIR = Path(__file__).resolve().parent


def mock_generate_chemshell_cluster(*args, **kwargs):
    # Copy the file to the output directory
    shutil.copytree(
        f"{FILE_DIR}/test_skzcam/ChemShell_Cluster.xyz_test",
        f"{FILE_DIR}/test_skzcam/ChemShell_Cluster.xyz",
    )


@pytest.fixture(autouse=True)
def patch_generate_chemshell_cluster(monkeypatch):
    monkeypatch.setattr(
        "quacc.atoms.skzcam.generate_chemshell_cluster", mock_generate_chemshell_cluster
    )
