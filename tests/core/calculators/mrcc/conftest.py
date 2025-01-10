from __future__ import annotations

from pathlib import Path

import pytest

FILE_DIR = Path(__file__).parent
MRCC_DIR = Path(FILE_DIR, "mrcc_run")


def mock_execute(self, directory, *args, **kwargs):
    import gzip

    with (
        gzip.open(MRCC_DIR / "mrcc.out.gz", "rb") as f,
        open(Path(directory, "mrcc.out"), "wb") as out,
    ):
        out.write(f.read())


@pytest.fixture(autouse=True)
def patch_execute(monkeypatch):
    from quacc.calculators.mrcc.mrcc import MrccTemplate

    monkeypatch.setattr(MrccTemplate, "execute", mock_execute)
