from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

FILE_DIR = Path(__file__).parent


def mock_run_chemshell(slab, slab_center_idx, atom_oxi_states, filepath, **kwargs):
    with (
        gzip.open(
            Path(FILE_DIR, "skzcam_files", "REF_ChemShell_Cluster.xyz.gz"), "rb"
        ) as f_in,
        Path(filepath, "ChemShell_Cluster.xyz").open(mode="wb") as f_out,
    ):
        shutil.copyfileobj(f_in, f_out)


@pytest.fixture(autouse=True)
def patch_run_chemshell(monkeypatch):
    from quacc.atoms.skzcam import CreateSKZCAMClusters

    monkeypatch.setattr(CreateSKZCAMClusters, "run_chemshell", mock_run_chemshell)
