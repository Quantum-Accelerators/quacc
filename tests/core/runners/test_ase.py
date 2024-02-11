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

from quacc import SETTINGS
from quacc.runners.ase import run_calc, run_opt, run_vib

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = True

DEFAULT_SETTINGS = SETTINGS.model_copy()


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
    SETTINGS.RESULTS_DIR = tmp_path

    prep_files()

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    new_atoms = run_calc(atoms, copy_files=["test_file.txt"])
    results_dir = _find_results_dir()

    assert atoms.calc.results is not None
    assert new_atoms.calc.results is not None
    assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
    assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))
    assert np.array_equal(new_atoms.get_positions(), atoms.get_positions()) is True
    assert np.array_equal(new_atoms.cell.array, atoms.cell.array) is True

    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR


def test_run_calc_no_gzip(tmp_path, monkeypatch):
    SETTINGS.RESULTS_DIR = tmp_path

    monkeypatch.chdir(tmp_path)
    prep_files()

    SETTINGS.GZIP_FILES = False

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    new_atoms = run_calc(atoms, copy_files=["test_file.txt"])
    results_dir = _find_results_dir()

    assert atoms.calc.results is not None
    assert new_atoms.calc.results is not None
    assert os.path.exists(os.path.join(results_dir, "test_file.txt"))
    assert not os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))
    assert np.array_equal(new_atoms.get_positions(), atoms.get_positions()) is True
    assert np.array_equal(new_atoms.cell.array, atoms.cell.array) is True
    SETTINGS.GZIP_FILES = DEFAULT_SETTINGS.GZIP_FILES
    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR


def test_run_opt1(tmp_path, monkeypatch):
    SETTINGS.RESULTS_DIR = tmp_path
    monkeypatch.chdir(tmp_path)
    prep_files()

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    dyn = run_opt(atoms, copy_files=["test_file.txt"])
    traj = dyn.traj_atoms
    results_dir = _find_results_dir()

    assert traj[-1].calc.results is not None
    assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
    assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))
    assert np.array_equal(traj[-1].get_positions(), atoms.get_positions()) is False
    assert np.array_equal(traj[-1].cell.array, atoms.cell.array) is True
    assert dyn.todict().get("restart")
    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR


def test_run_opt2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    dyn = run_opt(
        atoms,
        optimizer=BFGS,
        copy_files=["test_file.txt"],
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None

    dyn = run_opt(
        traj[-1],
        optimizer=BFGSLineSearch,
        copy_files=["test_file.txt"],
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None


def test_run_scipy_opt(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    dyn = run_opt(atoms, optimizer=SciPyFminBFGS)
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.todict().get("restart") is None


def test_run_vib(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prep_files()

    o2 = molecule("O2")
    o2.calc = LennardJones()
    vib = run_vib(o2, copy_files=["test_file.txt"])
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
        run_calc(atoms, copy_files=["test_file.txt"])

    # No file again
    with caplog.at_level(logging.WARNING):
        run_opt(atoms, copy_files=["test_file.txt"])

    # No trajectory kwarg
    with pytest.raises(ValueError):
        run_opt(
            atoms,
            optimizer=BFGSLineSearch,
            optimizer_kwargs={"restart": None, "trajectory": "test.traj"},
        )


def test_unique_workdir(tmp_path, monkeypatch):
    SETTINGS.CREATE_UNIQUE_DIR = True
    SETTINGS.RESULTS_DIR = tmp_path
    monkeypatch.chdir(tmp_path)
    prep_files()

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    run_calc(atoms, copy_files=["test_file.txt"])
    results_dir = _find_results_dir()
    assert atoms.calc.results is not None
    assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
    assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))

    SETTINGS.CREATE_UNIQUE_DIR = False
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    run_calc(atoms, copy_files=["test_file.txt"])
    results_dir = _find_results_dir()
    assert atoms.calc.results is not None
    assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
    assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))

    SETTINGS.CREATE_UNIQUE_DIR = DEFAULT_SETTINGS.CREATE_UNIQUE_DIR
    SETTINGS.RESULTS_DIR = DEFAULT_SETTINGS.RESULTS_DIR
