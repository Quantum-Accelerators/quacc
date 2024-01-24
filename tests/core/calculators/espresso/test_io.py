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
    assert (0, 0, 0) in results
    assert np.unique(results[(0, 0, 0)]["freqs"]).shape[0] == 1
    assert np.unique(results[(0, 0, 0)]["freqs"])[0] == 0.173268
    assert len(results[(0, 0, 0)]["eqpoints"]) == 1
    assert results[(0, 0, 0)]["atoms"].symbols == ["Al"]

    assert (0.75, -0.25, 0.75) in results
    assert np.unique(results[(0.75, -0.25, 0.75)]["freqs"]).shape[0] == 3
    assert np.unique(results[(0.75, -0.25, 0.75)]["freqs"])[2] == 8.791383
    assert len(results[(0.75, -0.25, 0.75)]["eqpoints"]) == 24
    assert results[(0.75, -0.25, 0.75)]["atoms"].symbols == ["Al"]


def test_dos_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    template = EspressoTemplate(binary="dos")
    template.outfiles = {"fildos": Path(RUN_PATH, "test.dos")}
    results = template.read_results(directory=".")
    assert results["test.dos"]["fermi"] == pytest.approx(5.98)
