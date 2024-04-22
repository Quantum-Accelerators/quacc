from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

FILE_DIR = Path(__file__).parent


def mock_generate_chemshell_cluster(*args, filepath, **kwargs):
    # Copy the file to the output directory
    with (
        gzip.open(
            Path(FILE_DIR, "skzcam_files", "REF_ChemShell_cluster.xyz.gz"), "rb"
        ) as f_in,
        Path(filepath, "ChemShell_cluster.xyz").open(mode="wb") as f_out,
    ):
        shutil.copyfileobj(f_in, f_out)


@pytest.fixture(autouse=True)
def patch_generate_chemshell_cluster(monkeypatch):
    from quacc.atoms import skzcam

    monkeypatch.setattr(
        skzcam, "generate_chemshell_cluster", mock_generate_chemshell_cluster
    )
