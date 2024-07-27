from __future__ import annotations

import gzip
import os

import pytest

from quacc.schemas.vasp import bader_runner, chargemol_runner


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
        "cm5": {"partial_charges": [1.0]},
    }


@pytest.fixture(autouse=True)
def patch_pop_analyses(monkeypatch):
    # Monkeypatch the Bader and Chargemol analyses so they don't run via pytest
    monkeypatch.setattr(
        "quacc.schemas.vasp.bader_analysis_from_path", mock_bader_analysis
    )
    monkeypatch.setattr("quacc.schemas.vasp.ChargemolAnalysis", mock_chargemol_analysis)


def prep_files():
    for f in ["CHGCAR", "POTCAR"]:
        with open(f, "w") as w:
            w.write("test")
    for f in ["AECCAR0.gz", "AECCAR2.gz"]:
        with gzip.open(f, "wb") as gf:
            gf.write(b"test")


def test_run_bader(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prep_files()

    bader_stats = bader_runner(tmp_path)
    assert bader_stats["min_dist"] == [1.0]
    assert bader_stats["partial_charges"] == [1.0]
    assert bader_stats["spin_moments"] == [0.0]
    assert bader_stats.get("charge") is None
    assert bader_stats.get("charge_transfer") is None
    assert bader_stats.get("reference_used") is None
    assert bader_stats.get("magmom") is None


def test_bader_erorr(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    with pytest.raises(FileNotFoundError):
        bader_runner(tmp_path)
    with open("CHGCAR", "w") as w:
        w.write("test")


def test_run_chargemol(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prep_files()

    chargemol_stats = chargemol_runner(path=tmp_path, atomic_densities_path=tmp_path)
    assert chargemol_stats["ddec"]["partial_charges"] == [1.0]
    assert chargemol_stats["ddec"]["spin_moments"] == [0.0]


def test_chargemol_erorr(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prep_files()

    with pytest.raises(
        OSError, match="DDEC6_ATOMIC_DENSITIES_DIR environment variable not defined"
    ):
        chargemol_runner(tmp_path)

    os.remove("CHGCAR")
    with pytest.raises(FileNotFoundError):
        chargemol_runner(tmp_path, atomic_densities_path=tmp_path)
    with open("CHGCAR", "w") as w:
        w.write("test")
