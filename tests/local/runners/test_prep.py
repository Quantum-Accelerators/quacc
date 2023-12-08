import os
from pathlib import Path

from quacc.runners.prep import calc_cleanup, calc_setup


def make_files():
    with open("file1.txt", "w") as f:
        f.write("file1")
    with open("file2.txt", "w") as f:
        f.write("file2")


def test_calc_setup(tmp_path, monkeypatch):
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.model_copy()
    monkeypatch.chdir(tmp_path)

    make_files()
    SETTINGS.SCRATCH_DIR = tmp_path

    tmpdir, results_dir = calc_setup()

    assert tmpdir.is_dir()
    assert "-tmp" in str(tmpdir)
    assert results_dir.is_dir()
    assert str(SETTINGS.RESULTS_DIR) in str(results_dir)
    if os.name != "nt":
        assert Path(SETTINGS.RESULTS_DIR, f"{tmpdir.name}-symlink").is_symlink()
    else:
        assert not Path(SETTINGS.RESULTS_DIR, f"{tmpdir.name}-symlink").is_symlink()
    assert "file1.txt" not in os.listdir(tmpdir)
    assert "file2.txt" not in os.listdir(tmpdir)
    assert Path.cwd() == tmpdir

    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR


def test_calc_setup_v2(tmp_path, monkeypatch):
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.model_copy()
    monkeypatch.chdir(tmp_path)

    make_files()
    SETTINGS.SCRATCH_DIR = None

    tmpdir, results_dir = calc_setup(copy_files=["file1.txt", "tmp_dir"])

    assert tmpdir.is_dir()
    assert "-tmp" in str(tmpdir)
    assert results_dir.is_dir()
    assert str(SETTINGS.RESULTS_DIR) in str(results_dir)
    assert not Path(SETTINGS.RESULTS_DIR, f"{tmpdir.name}-symlink").is_symlink()
    assert "file1.txt" in os.listdir(tmpdir)
    assert "file2.txt" not in os.listdir(tmpdir)
    assert Path.cwd() == tmpdir

    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR
