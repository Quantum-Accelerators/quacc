from __future__ import annotations

import os
from pathlib import Path
from shutil import rmtree

import pytest
from ase.build import bulk
from ase.calculators.emt import EMT

from quacc import change_settings
from quacc.runners.ase import Runner

jf = pytest.importorskip("jobflow")

CURRENT_DIR = Path(__file__).resolve().parent


def prep_files():
    # Make some test files to play with
    if not os.path.exists(CURRENT_DIR / "test_calc"):
        os.mkdir(CURRENT_DIR / "test_calc")
    with open(CURRENT_DIR / "test_file.txt", "w") as f:
        f.write("test")


def teardown_function():
    if os.path.exists(CURRENT_DIR / "test_calc"):
        rmtree(CURRENT_DIR / "test_calc", ignore_errors=True)
    for f in ["test_file.txt", "test_file.txt.gz"]:
        if os.path.exists(CURRENT_DIR / f):
            os.remove(CURRENT_DIR / f)


def test_run_calc(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    prep_files()

    @jf.job
    def test_job():
        with change_settings({"RESULTS_DIR": tmp_path, "GZIP_FILES": False}):
            atoms = bulk("Cu") * (2, 1, 1)
            atoms[0].position += 0.1
            assert "test_file.txt" not in os.listdir(tmp_path)
            Runner(atoms, EMT(), copy_files={CURRENT_DIR: "test_file.txt"}).run_calc()
            assert "test_file.txt" in os.listdir(tmp_path)

    jf.run_locally(test_job(), ensure_success=True, create_folders=False)


def test_run_calc2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    prep_files()

    @jf.job
    def test_job():
        with change_settings(
            {"RESULTS_DIR": tmp_path, "GZIP_FILES": False, "CREATE_UNIQUE_DIR": False}
        ):
            atoms = bulk("Cu") * (2, 1, 1)
            atoms[0].position += 0.1
            assert "test_file.txt" not in os.listdir(tmp_path)
            Runner(atoms, EMT(), copy_files={CURRENT_DIR: "test_file.txt"}).run_calc()
            assert "test_file.txt" in os.listdir(tmp_path)

    jf.run_locally(test_job(), ensure_success=True, create_folders=True)
