from __future__ import annotations

import os
from pathlib import Path

import pytest
from ase.build import bulk
from ase.calculators.emt import EMT

jf = pytest.importorskip("jobflow")


from quacc import change_settings, job
from quacc.runners.prep import calc_setup

CURRENT_DIR = Path(__file__).resolve().parent


def make_files():
    with open(CURRENT_DIR / "file1.txt", "w") as f:
        f.write("file1")
    with open(CURRENT_DIR / "file2.txt", "w") as f:
        f.write("file2")


def teardown_function():
    for file in ["file1.txt", "file2.txt"]:
        with contextlib.suppress(FileNotFoundError):
            os.remove(CURRENT_DIR / file)


@pytest.mark.parametrize(
    "copy_files",
    [{CURRENT_DIR: ["file1.txt"]}, {CURRENT_DIR: "file1.txt"}, {CURRENT_DIR: "file1*"}],
)
def test_calc_setup_v1(tmp_path, monkeypatch, copy_files):
    monkeypatch.chdir(tmp_path)
    make_files()

    @job
    def test_job():
        with change_settings({"CREATE_UNIQUE_DIR": True}):
            atoms = bulk("Cu")
            atoms.calc = EMT()

            tmpdir, _ = calc_setup(atoms, copy_files=copy_files)

            assert tmpdir.is_dir()
            assert "tmp" in str(tmpdir)
            assert "file1.txt" in os.listdir(tmpdir)
            assert "file2.txt" not in os.listdir(tmpdir)

    jf.run_locally(test_job(), ensure_success=True, create_folders=False)


@pytest.mark.parametrize(
    "copy_files",
    [{CURRENT_DIR: ["file1.txt"]}, {CURRENT_DIR: "file1.txt"}, {CURRENT_DIR: "file1*"}],
)
def test_calc_setup_v2(tmp_path, monkeypatch, copy_files):
    monkeypatch.chdir(tmp_path)
    make_files()

    @job
    def test_job():
        with change_settings({"CREATE_UNIQUE_DIR": False}):
            atoms = bulk("Cu")
            atoms.calc = EMT()

            tmpdir, _ = calc_setup(atoms, copy_files=copy_files)

            assert tmpdir.is_dir()
            assert "file1.txt" in os.listdir(tmpdir)
            assert "file2.txt" not in os.listdir(tmpdir)

    jf.run_locally(test_job(), ensure_success=True, create_folders=False)


@pytest.mark.parametrize(
    "copy_files",
    [{CURRENT_DIR: ["file1.txt"]}, {CURRENT_DIR: "file1.txt"}, {CURRENT_DIR: "file1*"}],
)
def test_calc_setup_v3(tmp_path, monkeypatch, copy_files):
    monkeypatch.chdir(tmp_path)
    make_files()

    @job
    def test_job():
        with change_settings({"CREATE_UNIQUE_DIR": False}):
            atoms = bulk("Cu")
            atoms.calc = EMT()
            tmpdir, _ = calc_setup(atoms, copy_files=copy_files)

            assert tmpdir.is_dir()
            assert "file1.txt" in os.listdir(tmpdir)
            assert "file2.txt" not in os.listdir(tmpdir)

    jf.run_locally(test_job(), ensure_success=True, create_folders=True)
