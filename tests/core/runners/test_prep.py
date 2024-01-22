import os
from pathlib import Path

from ase.build import bulk
from ase.calculators.emt import EMT

from quacc.runners.prep import calc_cleanup, calc_setup


def make_files():
    with open("file1.txt", "w") as f:
        f.write("file1")
    with open("file2.txt", "w") as f:
        f.write("file2")


def make_files2():
    p = Path("quacc-tmp-1234")
    os.makedirs(p)


def test_calc_setup(tmp_path, monkeypatch):
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.model_copy()
    monkeypatch.chdir(tmp_path)

    make_files()
    SETTINGS.SCRATCH_DIR = tmp_path
    atoms = bulk("Cu")
    atoms.calc = EMT()

    tmpdir, results_dir = calc_setup(atoms)

    assert tmpdir.is_dir()
    assert "tmp" in str(tmpdir)
    assert results_dir.name == tmpdir.name.split("tmp-")[-1]
    assert str(SETTINGS.RESULTS_DIR) in str(results_dir)
    if os.name != "nt":
        assert Path(SETTINGS.RESULTS_DIR, f"symlink-{tmpdir.name}").is_symlink()
    else:
        assert not Path(SETTINGS.RESULTS_DIR, f"symlink-{tmpdir.name}").exists()
    assert "file1.txt" not in os.listdir(tmpdir)
    assert "file2.txt" not in os.listdir(tmpdir)
    assert Path.cwd() == tmpdir
    assert Path(atoms.calc.directory) == tmpdir

    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR


def test_calc_setup_v2(tmp_path, monkeypatch):
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.model_copy()
    monkeypatch.chdir(tmp_path)

    make_files()
    SETTINGS.SCRATCH_DIR = None
    atoms = bulk("Cu")
    atoms.calc = EMT()

    tmpdir, results_dir = calc_setup(atoms, copy_files=["file1.txt", "tmp_dir"])

    assert tmpdir.is_dir()
    assert "tmp" in str(tmpdir)
    assert results_dir.name == tmpdir.name.split("tmp-")[-1]
    assert str(SETTINGS.RESULTS_DIR) in str(results_dir)
    assert not Path(SETTINGS.RESULTS_DIR, f"symlink-{tmpdir.name}").exists()
    assert "file1.txt" in os.listdir(tmpdir)
    assert "file2.txt" not in os.listdir(tmpdir)
    assert Path.cwd() == tmpdir

    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR


def test_calc_cleanup(tmp_path, monkeypatch):
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.model_copy()
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    atoms.calc = EMT()

    make_files2()
    SETTINGS.SCRATCH_DIR = tmp_path

    p = Path(Path.cwd(), "quacc-tmp-1234").resolve()
    assert p.is_dir()
    calc_cleanup(atoms, p, SETTINGS.RESULTS_DIR)
    assert not p.exists()
    assert Path.cwd() == SETTINGS.RESULTS_DIR
    assert Path(atoms.calc.directory) == SETTINGS.RESULTS_DIR

    SETTINGS.SCRATCH_DIR = DEFAULT_SETTINGS.SCRATCH_DIR
