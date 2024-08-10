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
from ase.io import read
from ase.optimize import BFGS, BFGSLineSearch
from ase.optimize.sciopt import SciPyFminBFGS

from quacc import JobFailure, change_settings, get_settings
from quacc.runners._base import BaseRunner
from quacc.runners.ase import Runner

LOGGER = logging.getLogger(__name__)
LOGGER.propagate = True


def _find_results_dir():
    search_dir = get_settings().RESULTS_DIR
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
    results_dir = get_settings().RESULTS_DIR
    if os.path.exists(os.path.join(results_dir, "test_calc")):
        rmtree(os.path.join(results_dir, "test_calc"), ignore_errors=True)
    for f in ["test_file.txt", "test_file.txt.gz"]:
        if os.path.exists(os.path.join(results_dir, f)):
            os.remove(os.path.join(results_dir, f))


def test_base_runner(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    atoms.calc = EMT()
    br = BaseRunner(atoms)

    br.setup()
    assert "tmp" in str(br.tmpdir)
    assert br.tmpdir.exists()
    assert "tmp" not in str(br.job_results_dir)
    assert not br.job_results_dir.exists()
    assert Path(br.atoms.calc.directory) == br.tmpdir

    br.cleanup()
    assert not br.tmpdir.exists()
    assert br.job_results_dir.exists()
    assert Path(br.atoms.calc.directory) == br.job_results_dir


def test_base_runner2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    br = BaseRunner()

    br.setup()
    assert "tmp" in str(br.tmpdir)
    assert br.tmpdir.exists()
    assert "tmp" not in str(br.job_results_dir)
    assert not br.job_results_dir.exists()

    br.cleanup()
    assert not br.tmpdir.exists()
    assert br.job_results_dir.exists()


def test_run_calc(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    prep_files()

    with change_settings({"RESULTS_DIR": tmp_path}):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.1

        new_atoms = Runner(
            atoms, EMT(), copy_files={Path(): "test_file.txt"}
        ).run_calc()
        results_dir = _find_results_dir()

        assert atoms.calc is None
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

        new_atoms = Runner(
            atoms, EMT(), copy_files={Path(): "test_file.txt"}
        ).run_calc()
        results_dir = _find_results_dir()

        assert atoms.calc is None
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

        dyn = Runner(atoms, EMT(), copy_files={Path(): "test_file.txt"}).run_opt()
        traj = read(dyn.trajectory.filename, index=":")
        results_dir = _find_results_dir()

        assert atoms.calc is None
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

    dyn = Runner(atoms, EMT(), copy_files={Path(): "test_file.txt"}).run_opt(
        optimizer=BFGS, optimizer_kwargs={"restart": None}
    )
    traj = read(dyn.trajectory.filename, index=":")
    assert traj[-1].calc.results is not None

    dyn = Runner(traj[-1], EMT(), copy_files={Path(): "test_file.txt"}).run_opt(
        optimizer=BFGSLineSearch, optimizer_kwargs={"restart": None}
    )
    traj = read(dyn.trajectory.filename, index=":")
    assert traj[-1].calc.results is not None


def test_run_scipy_opt(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1

    dyn = Runner(atoms, EMT()).run_opt(optimizer=SciPyFminBFGS)
    traj = read(dyn.trajectory.filename, index=":")
    assert traj[-1].calc.results is not None
    assert dyn.todict().get("restart") is None


def test_run_vib(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prep_files()

    atoms = molecule("O2")
    vib = Runner(atoms, LennardJones(), copy_files={Path(): "test_file.txt"}).run_vib()
    results_dir = _find_results_dir()

    assert atoms.calc is None
    assert np.real(vib.get_frequencies()[-1]) == pytest.approx(255.6863883406967)
    assert np.array_equal(vib.atoms.get_positions(), atoms.get_positions()) is True
    assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
    assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))


def test_bad_runs(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    # No file
    with caplog.at_level(logging.WARNING):
        Runner(atoms, EMT(), copy_files={Path(): "test_file.txt"}).run_calc()
    assert "Cannot find file" in caplog.text

    # No file again
    with caplog.at_level(logging.WARNING):
        Runner(atoms, EMT(), copy_files={Path(): "test_file.txt"}).run_opt()
    assert "Cannot find file" in caplog.text

    # No trajectory kwarg
    with pytest.raises(
        ValueError, match="Quacc does not support setting the `trajectory` kwarg"
    ):
        Runner(atoms, EMT()).run_opt(
            optimizer=BFGSLineSearch,
            optimizer_kwargs={"restart": None, "trajectory": "test.traj"},
        )


def test_unique_workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    prep_files()

    with change_settings({"CREATE_UNIQUE_DIR": True, "RESULTS_DIR": tmp_path}):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.1

        Runner(atoms, EMT(), copy_files={Path(): "test_file.txt"}).run_calc()
        results_dir = _find_results_dir()
        assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
        assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))

    with change_settings({"CREATE_UNIQUE_DIR": False, "RESULTS_DIR": tmp_path}):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.1

        Runner(atoms, EMT(), copy_files={Path(): "test_file.txt"}).run_calc()
        results_dir = _find_results_dir()
        assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
        assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))


def test_fn_hook(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def fn_hook(dyn):
        if dyn.atoms:
            raise ValueError("Test error")

    with pytest.raises(JobFailure, match="Calculation failed!") as err:
        Runner(bulk("Cu"), EMT()).run_opt(fn_hook=fn_hook)
    with pytest.raises(RuntimeError, match="Test error"):
        raise err.value.parent_error
