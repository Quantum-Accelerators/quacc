from __future__ import annotations

import os
from pathlib import Path

import pytest
from ase.build import bulk
from ase.calculators.emt import EMT

from quacc import JobFailure, change_settings, get_settings
from quacc.runners.prep import calc_cleanup, calc_setup, terminate


def make_files():
    with open("file1.txt", "w") as f:
        f.write("file1")
    with open("file2.txt", "w") as f:
        f.write("file2")


def make_files2():
    p = Path("quacc-tmp-1234")
    os.makedirs(p)


def make_files3():
    Path("saved").mkdir()
    with open("saved/file1.txt", "w") as f:
        f.write("file1")
    with open("saved/file2.txt", "w") as f:
        f.write("file2")


def test_calc_setup(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    make_files()

    with change_settings({"SCRATCH_DIR": tmp_path}):
        atoms = bulk("Cu")
        atoms.calc = EMT()
        settings = get_settings()

        tmpdir, results_dir = calc_setup(atoms)

        assert tmpdir.is_dir()
        assert "tmp" in str(tmpdir)
        assert results_dir.name == tmpdir.name.split("tmp-")[-1]
        assert str(settings.RESULTS_DIR) in str(results_dir)
        if os.name != "nt":
            assert Path(settings.RESULTS_DIR, f"symlink-{tmpdir.name}").is_symlink()
        else:
            assert not Path(settings.RESULTS_DIR, f"symlink-{tmpdir.name}").exists()
        assert "file1.txt" not in os.listdir(tmpdir)
        assert "file2.txt" not in os.listdir(tmpdir)
        assert Path(atoms.calc.directory) == tmpdir


@pytest.mark.parametrize(
    "copy_files", [{Path(): ["file1.txt"]}, {Path(): "file1.txt"}, {Path(): "file1*"}]
)
def test_calc_setup_v2(tmp_path, monkeypatch, copy_files):
    monkeypatch.chdir(tmp_path)
    make_files()

    with change_settings({"SCRATCH_DIR": None}):
        atoms = bulk("Cu")
        atoms.calc = EMT()
        settings = get_settings()

        tmpdir, results_dir = calc_setup(atoms, copy_files=copy_files)

        assert tmpdir.is_dir()
        assert "tmp" in str(tmpdir)
        assert results_dir.name == tmpdir.name.split("tmp-")[-1]
        assert str(settings.RESULTS_DIR) in str(results_dir)
        assert not Path(settings.RESULTS_DIR, f"symlink-{tmpdir.name}").exists()
        assert "file1.txt" in os.listdir(tmpdir)
        assert "file2.txt" not in os.listdir(tmpdir)


@pytest.mark.parametrize(
    "copy_files",
    [
        {Path("saved"): "file1.txt"},
        {"saved": "file1.txt"},
        {"saved": "file1*"},
        {"saved": ["file1.txt"]},
        {"saved": ["file1*"]},
    ],
)
def test_calc_setup_v3(tmp_path, monkeypatch, copy_files):
    monkeypatch.chdir(tmp_path)
    make_files3()

    with change_settings({"SCRATCH_DIR": None}):
        atoms = bulk("Cu")
        atoms.calc = EMT()
        settings = get_settings()

        tmpdir, results_dir = calc_setup(atoms, copy_files=copy_files)

        assert tmpdir.is_dir()
        assert "tmp" in str(tmpdir)
        assert results_dir.name == tmpdir.name.split("tmp-")[-1]
        assert str(settings.RESULTS_DIR) in str(results_dir)
        assert not Path(settings.RESULTS_DIR, f"symlink-{tmpdir.name}").exists()
        assert "saved" not in os.listdir(tmpdir)
        assert "file1.txt" in os.listdir(tmpdir)
        assert "file2.txt" not in os.listdir(tmpdir)


@pytest.mark.parametrize("copy_files", ["saved", Path("saved")])
def test_calc_setup_v3_2(tmp_path, monkeypatch, copy_files):
    monkeypatch.chdir(tmp_path)
    make_files3()

    with change_settings({"SCRATCH_DIR": None}):
        atoms = bulk("Cu")
        atoms.calc = EMT()
        settings = get_settings()

        tmpdir, results_dir = calc_setup(atoms, copy_files=copy_files)

        assert tmpdir.is_dir()
        assert "tmp" in str(tmpdir)
        assert results_dir.name == tmpdir.name.split("tmp-")[-1]
        assert str(settings.RESULTS_DIR) in str(results_dir)
        assert not Path(settings.RESULTS_DIR, f"symlink-{tmpdir.name}").exists()
        assert "saved" not in os.listdir(tmpdir)
        assert "file1.txt" in os.listdir(tmpdir)
        assert "file2.txt" in os.listdir(tmpdir)


@pytest.mark.parametrize("copy_files", [{"saved": "*"}])
def test_calc_setup_v4(tmp_path, monkeypatch, copy_files):
    monkeypatch.chdir(tmp_path)
    make_files3()

    with change_settings({"SCRATCH_DIR": None}):
        atoms = bulk("Cu")
        atoms.calc = EMT()
        settings = get_settings()

        tmpdir, results_dir = calc_setup(atoms, copy_files=copy_files)

        assert tmpdir.is_dir()
        assert "tmp" in str(tmpdir)
        assert results_dir.name == tmpdir.name.split("tmp-")[-1]
        assert str(settings.RESULTS_DIR) in str(results_dir)
        assert not Path(settings.RESULTS_DIR, f"symlink-{tmpdir.name}").exists()
        assert "file1.txt" in os.listdir(tmpdir)
        assert "file2.txt" in os.listdir(tmpdir)
        assert "saved" not in os.listdir(tmpdir)


def test_calc_cleanup(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    atoms.calc = EMT()

    make_files2()

    with change_settings({"SCRATCH_DIR": tmp_path}):
        settings = get_settings()
        p = Path(tmp_path, "quacc-tmp-1234").resolve()
        assert p.is_dir()
        calc_cleanup(atoms, p, settings.RESULTS_DIR)
        assert not p.exists()
        assert Path(atoms.calc.directory) == settings.RESULTS_DIR

    with pytest.raises(
        ValueError, match="quacc does not appear to be a tmpdir"
    ):
        calc_cleanup(atoms, "quacc", get_settings().RESULTS_DIR)


def test_terminate(tmp_path):
    p = tmp_path / "tmp-quacc-1234"
    os.mkdir(p)
    with pytest.raises(JobFailure, match="Calculation failed!") as err:
        terminate(p, ValueError("moo"))
    with pytest.raises(ValueError, match="moo"):
        raise err.value.parent_error
    assert err.value.directory == tmp_path / "failed-quacc-1234"
    assert not p.exists()
    assert Path(tmp_path, "failed-quacc-1234").exists()
