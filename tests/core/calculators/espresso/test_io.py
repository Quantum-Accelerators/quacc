import gzip
from io import StringIO
from pathlib import Path

import numpy as np
import pytest
from ase.io.espresso import read_espresso_ph

from quacc.calculators.espresso.espresso import EspressoTemplate

RUN_PATH = Path(__file__).parent / "test_runs"


def _read_gzip_file(file_path):
    with gzip.open(file_path, "rt", encoding="utf-8") as file:
        content = file.read()
    return content


@pytest.mark.parametrize("i", range(1, 66, 5))
def test_read_espresso_ph_all(i):
    file_contents = _read_gzip_file(RUN_PATH / f"run{i}.txt.gz")
    fd = StringIO(file_contents)
    read_espresso_ph(fd)


def test_read_espresso_ph_1():
    file_contents = _read_gzip_file(RUN_PATH / "run1.txt.gz")
    fd = StringIO(file_contents)
    results = read_espresso_ph(fd)

    assert len(results) == 8
    assert np.unique(results[1]["freqs"]).shape[0] == 1
    assert np.unique(results[1]["freqs"])[0] == 0.173268
    assert len(results[1]["eqpoints"]) == 1
    assert results[1]["atoms"].symbols == ["Al"]

    assert np.unique(results[5]["freqs"]).shape[0] == 3
    assert np.unique(results[5]["freqs"])[2] == 8.791383
    assert len(results[5]["eqpoints"]) == 24
    assert results[5]["atoms"].symbols == ["Al"]


def test_dos_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    template = EspressoTemplate(binary="dos")
    template.outfiles = {"fildos": Path(RUN_PATH, "test.dos")}
    results = template.read_results(directory=".")
    assert results["test_dos"]["fermi"] == pytest.approx(5.98)
