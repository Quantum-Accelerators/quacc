import gzip
import os
from shutil import rmtree

import pytest

from quacc.util.pop_analysis import bader_runner, chargemol_runner


def mock_bader_analysis(*args, **kwargs):
    return {
        "min_dist": [1.0],
        "charge": [1.0],
        "atomic_volume": [1.0],
        "vacuum_charge": 1.0,
        "vacuum_volume": 1.0,
        "bader_version": 1.0,
        "reference_used": [0.0],
        "charge_transfer": [-1.0],
        "magmom": [0.0],
    }


def mock_chargemol_analysis(*args, **kwargs):
    return {
        "ddec": {
            "partial_charges": [1.0],
            "spin_moments": [0.0],
            "dipoles": [1.0],
            "rsquared_moments": [1.0],
            "rcubed_moments": [1.0],
            "rfourth_moments": [1.0],
        },
        "cm5": {
            "partial_charges": [1.0],
        },
    }


@pytest.fixture(autouse=True)
def patch_pop_analyses(monkeypatch):
    # Monkeypatch the Bader and Chargemol analyses so they don't run via pytest
    monkeypatch.setattr(
        "quacc.util.pop_analysis.bader_analysis_from_path",
        mock_bader_analysis,
    )
    monkeypatch.setattr(
        "quacc.util.pop_analysis.ChargemolAnalysis",
        mock_chargemol_analysis,
    )


CWD = os.getcwd()


def setup_module():
    if not os.path.exists("rundir"):
        os.mkdir("rundir")
    os.chdir("rundir")
    for f in ["CHGCAR", "POTCAR"]:
        with open(f, "w") as w:
            w.write("test")
    for f in ["AECCAR0.gz", "AECCAR2.gz"]:
        with gzip.open(f, "wb") as gf:
            gf.write("test".encode())


def teardown_module():
    for f in ["CHGCAR", "AECCAR0.gz", "AECCAR2.gz", "POTCAR"]:
        if os.path.isfile(f):
            os.remove(f)
    if os.path.exists("rundir"):
        rmtree("rundir")
    os.chdir(CWD)


def test_run_bader():
    bader_stats = bader_runner()
    assert bader_stats["min_dist"] == [1.0]
    assert bader_stats["partial_charges"] == [1.0]
    assert bader_stats["spin_moments"] == [0.0]
    assert bader_stats.get("charge") is None
    assert bader_stats.get("charge_transfer") is None
    assert bader_stats.get("reference_used") is None
    assert bader_stats.get("magmom") is None


def test_bader_erorr():
    os.remove("CHGCAR")
    with pytest.raises(FileNotFoundError):
        bader_runner()
    with open("CHGCAR", "w") as w:
        w.write("test")


def test_run_chargemol():
    chargemol_stats = chargemol_runner(path=".", atomic_densities_path=".")
    assert chargemol_stats["ddec"]["partial_charges"] == [1.0]
    assert chargemol_stats["ddec"]["spin_moments"] == [0.0]


def test_chargemol_erorr():
    with pytest.raises(ValueError):
        chargemol_runner()

    os.remove("CHGCAR")
    with pytest.raises(FileNotFoundError):
        chargemol_runner(atomic_densities_path=".")
    with open("CHGCAR", "w") as w:
        w.write("test")
