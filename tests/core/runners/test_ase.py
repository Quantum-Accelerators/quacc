from __future__ import annotations

import glob
import logging
import os
from pathlib import Path
from shutil import rmtree

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS, BFGSLineSearch
from ase.optimize.sciopt import SciPyFminBFGS

from quacc import SETTINGS, change_settings
from quacc import SETTINGS
from quacc.runners.ase import Runner

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = True


def _find_results_dir():
    search_dir = SETTINGS.RESULTS_DIR
    pattern = str(Path(search_dir, "quacc-*"))
    matching_dirs = glob.glob(pattern)
    most_recent_directory = max(matching_dirs, key=os.path.getmtime, default=None)
    return Path.cwd() if most_recent_directory is None else most_recent_directory


def prep_files():
    # Make some test files to play with
    if not os.path.exists("test_calc"):
        os.mkdir("test_calc")
    with open("test_file.txt", "w") as f:
        f.write("test")


def teardown_function():
    if os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, "test_calc")):
        rmtree(os.path.join(SETTINGS.RESULTS_DIR, "test_calc"), ignore_errors=True)
    for f in ["test_file.txt", "test_file.txt.gz"]:
        if os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, f)):
            os.remove(os.path.join(SETTINGS.RESULTS_DIR, f))


def test_run_calc(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    prep_files()

    with change_settings({"RESULTS_DIR": tmp_path}):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.1
        atoms.calc = EMT()

        new_atoms = Runner(atoms, copy_files={Path(): "test_file.txt"}).run_calc()
        results_dir = _find_results_dir()

        assert atoms.calc.results is not None
        assert new_atoms.calc.results is not None
        assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
        assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))
        assert np.array_equal(new_atoms.get_positions(), atoms.get_positions()) is True
        assert np.array_equal(new_atoms.cell.array, atoms.cell.array) is True


def test_run_calc_no_gzip(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prep_files()

    with change_settings({"RESULTS_DIR": tmp_path, "GZIP_FILES": False}):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.1
        atoms.calc = EMT()

        new_atoms = Runner(atoms, copy_files={Path(): "test_file.txt"}).run_calc()
        results_dir = _find_results_dir()

        assert atoms.calc.results is not None
        assert new_atoms.calc.results is not None
        assert os.path.exists(os.path.join(results_dir, "test_file.txt"))
        assert not os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))
        assert np.array_equal(new_atoms.get_positions(), atoms.get_positions()) is True
        assert np.array_equal(new_atoms.cell.array, atoms.cell.array) is True


def test_run_opt1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prep_files()

    with change_settings({"RESULTS_DIR": tmp_path}):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.1
        atoms.calc = EMT()

        dyn = Runner(atoms, copy_files={Path(): "test_file.txt"}).run_opt()
        traj = dyn.traj_atoms
        results_dir = _find_results_dir()

        assert traj[-1].calc.results is not None
        assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
        assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))
        assert np.array_equal(traj[-1].get_positions(), atoms.get_positions()) is False
        assert np.array_equal(traj[-1].cell.array, atoms.cell.array) is True
        assert dyn.todict().get("restart")


def test_run_opt2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    dyn = Runner(atoms, copy_files={Path(): "test_file.txt"}).run_opt(
        optimizer=BFGS,
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None

    dyn = Runner(traj[-1], copy_files={Path(): "test_file.txt"}).run_opt(
        optimizer=BFGSLineSearch,
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None


def test_run_scipy_opt(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    dyn = Runner(atoms).run_opt(optimizer=SciPyFminBFGS)
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.todict().get("restart") is None


def test_run_vib(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prep_files()

    o2 = molecule("O2")
    o2.calc = LennardJones()
    vib = Runner(o2, copy_files={Path(): "test_file.txt"}).run_vib()
    results_dir = _find_results_dir()

    assert np.real(vib.get_frequencies()[-1]) == pytest.approx(255.6863883406967)
    assert np.array_equal(vib.atoms.get_positions(), o2.get_positions()) is True
    assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
    assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))


def test_bad_runs(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    atoms.calc = EMT()

    # No file
    with caplog.at_level(logging.WARNING):
        Runner(atoms, copy_files={Path(): "test_file.txt"}).run_calc()
    assert "Cannot find file" in caplog.text

    # No file again
    with caplog.at_level(logging.WARNING):
        Runner(atoms, copy_files={Path(): "test_file.txt"}).run_opt()
    assert "Cannot find file" in caplog.text

    # No trajectory kwarg
    with pytest.raises(ValueError):
        Runner(atoms).run_opt(
            optimizer=BFGSLineSearch,
            optimizer_kwargs={"restart": None, "trajectory": "test.traj"},
        )


def test_unique_workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prep_files()

    with change_settings({"CREATE_UNIQUE_DIR": True, "RESULTS_DIR": tmp_path}):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.1
        atoms.calc = EMT()

        Runner(atoms, copy_files={Path(): "test_file.txt"}).run_calc()
        results_dir = _find_results_dir()
        assert atoms.calc.results is not None
        assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
        assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))

    with change_settings({"CREATE_UNIQUE_DIR": False, "RESULTS_DIR": tmp_path}):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.1
        atoms.calc = EMT()

        Runner(atoms, copy_files={Path(): "test_file.txt"}).run_calc()
        results_dir = _find_results_dir()
        assert atoms.calc.results is not None
        assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
        assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))


def test_fn_hook(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def fn_hook(dyn):
        if dyn.atoms:
            raise ValueError("Test error")

    atoms = bulk("Cu")
    atoms.calc = EMT()

    with pytest.raises(ValueError, match="Test error"):
        run_opt(atoms, fn_hook=fn_hook)
