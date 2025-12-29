from __future__ import annotations

import glob
import os
from importlib.util import find_spec
from logging import WARNING, getLogger
from pathlib import Path
from shutil import rmtree

import numpy as np
import pytest


@pytest.fixture(scope="module", autouse=True)
def set_seed():
    np.random.seed(42)  # noqa: NPY002


from ase import Atoms
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase.io import read
from ase.mep.neb import NEBOptimizer
from ase.optimize import BFGS, BFGSLineSearch
from ase.optimize.sciopt import SciPyFminBFGS

from quacc import JobFailure, change_settings, get_settings
from quacc.runners._base import BaseRunner
from quacc.runners.ase import Runner

has_geodesic_interpolate = bool(find_spec("geodesic_interpolate"))
test_files_path = Path(__file__).parent / "test_files"

LOGGER = getLogger(__name__)
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


@pytest.fixture
def setup_test_environment(tmp_path):
    reactant = Atoms(
        symbols="CCHHCHH",
        positions=[
            [1.4835950817281542, -1.0145410211301968, -0.13209027203235943],
            [0.8409564131524673, 0.018549610257914483, -0.07338809662321308],
            [-0.6399757891931867, 0.01763740851518944, 0.0581573443268891],
            [-1.0005576455546672, 1.0430257532387608, 0.22197240310602892],
            [1.402180736662139, 0.944112416574632, -0.12179540364365492],
            [-1.1216961389434357, -0.3883639833876232, -0.8769102842015071],
            [-0.9645026578514683, -0.6204201840686793, 0.9240543090678239],
        ],
    )

    product = Atoms(
        symbols="CCHHCHH",
        positions=[
            [1.348003553501624, 0.4819311116778978, 0.2752537177143993],
            [0.2386618286631742, -0.3433222966734429, 0.37705518940917926],
            [-0.9741307940518336, 0.07686022294949588, 0.08710778043683955],
            [-1.8314843503320921, -0.5547344604780035, 0.1639037492534953],
            [0.3801391040059668, -1.3793340533058087, 0.71035902765307],
            [1.9296265384257907, 0.622088341468767, 1.0901733942191298],
            [-1.090815880212625, 1.0965111343610956, -0.23791518420660265],
        ],
    )
    return reactant, product


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
    with caplog.at_level(WARNING):
        Runner(atoms, EMT(), copy_files={Path(): "test_file.txt"}).run_calc()
    assert "Cannot find file" in caplog.text
    caplog.clear()

    # No file again
    with caplog.at_level(WARNING):
        Runner(atoms, EMT(), copy_files={Path(): "test_file.txt"}).run_opt()
    assert "Cannot find file" in caplog.text
    caplog.clear()

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
    with pytest.raises(ValueError, match="Test error"):
        raise err.value.parent_error


def test_run_neb(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    geodesic_path = test_files_path / "geodesic_path.xyz"
    images = read(geodesic_path, index=":")

    neb_kwargs = {"method": "aseneb", "precon": None}
    dyn = Runner(images, EMT()).run_neb(optimizer=NEBOptimizer, neb_kwargs=neb_kwargs)
    traj = read(dyn.trajectory.filename, index=":")

    assert traj[-1].calc.results is not None
    assert not os.path.exists(tmp_path / "opt.log")


def test_run_neb2():
    geodesic_path = test_files_path / "geodesic_path.xyz"

    images = read(geodesic_path, index=":")

    with pytest.raises(
        ValueError, match="BFGSLineSearch is not allowed as optimizer with NEB"
    ):
        Runner(images, EMT()).run_neb(
            optimizer=BFGSLineSearch, neb_kwargs={"method": "aseneb", "precon": None}
        )


def test_run_neb_raises_value_error_for_trajectory_kwarg():
    images = [Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])]

    with pytest.raises(
        ValueError, match="Quacc does not support setting the `trajectory` kwarg"
    ):
        Runner(images, EMT()).run_neb(
            optimizer=NEBOptimizer, optimizer_kwargs={"trajectory": "some_traj.traj"}
        )
