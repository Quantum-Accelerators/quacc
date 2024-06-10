from __future__ import annotations

import glob
import logging
import os
from pathlib import Path
from shutil import rmtree

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase.mep.neb import NEBOptimizer
from ase.optimize import BFGS, BFGSLineSearch
from ase.optimize.sciopt import SciPyFminBFGS
from sella import Sella

from quacc import SETTINGS, change_settings, strip_decorator
from quacc.recipes.newtonnet.ts import relax_job
from quacc.runners.ase import (
    _geodesic_interpolate_wrapper,
    run_calc,
    run_opt,
    run_path_opt,
    run_vib,
)
from quacc.schemas.ase import summarize_path_opt_run

has_newtonnet = bool(find_spec("newtonnet"))
if has_newtonnet:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet

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


@pytest.fixture()
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
    current_file_path = Path(__file__).parent
    conf_path = (
        current_file_path / "../../../tests/core/recipes/newtonnet_recipes"
    ).resolve()
    NEWTONNET_CONFIG_PATH = conf_path / "config0.yml"
    NEWTONNET_MODEL_PATH = conf_path / "best_model_state.tar"
    SETTINGS.CHECK_CONVERGENCE = False
    calc_defaults = {
        "model_path": NEWTONNET_MODEL_PATH,
        "settings_path": NEWTONNET_CONFIG_PATH,
    }
    SETTINGS.NEWTONNET_MODEL_PATH = NEWTONNET_MODEL_PATH
    SETTINGS.NEWTONNET_CONFIG_PATH = NEWTONNET_CONFIG_PATH
    return reactant, product, calc_defaults


@pytest.mark.parametrize(
    (
        "method",
        "optimizer_class",
        "precon",
        "n_intermediate",
        "k",
        "max_steps",
        "fmax",
        "expected_logfile",
        "first_image_positions",
        "first_image_pot_energy",
        "first_image_forces",
        "second_images_positions",
        "index_ts",
        "pot_energy_ts",
        "forces_ts",
        "last_images_positions",
    ),
    [
        # ("aseneb", NEBOptimizer, None, 10, 0.1, 3, 1e-3, None),
        # ("aseneb", SciPyFminBFGS, None, 1000, 0.1, 3, 1e-3, "some_logdir",
        #   0.78503956131, -24.9895786292, -0.0017252843, 0.78017739462, 9, -19.946616164,
        #   -0.19927549, 0.51475535802),
        (
            "aseneb",
            NEBOptimizer,
            None,
            10,
            0.1,
            3,
            1e-3,
            "some_logdir",
            0.78503956131,
            -24.9895786292,
            -0.0017252843,
            0.78017739462,
            9,
            -19.946616164,
            -0.19927549,
            0.51475535802,
        )
    ],
)
def test_run_neb_method(
    setup_test_environment,
    tmp_path,
    method,
    optimizer_class,
    precon,
    n_intermediate,
    k,
    max_steps,
    fmax,
    expected_logfile,
    first_image_positions,
    first_image_pot_energy,
    first_image_forces,
    second_images_positions,
    index_ts,
    pot_energy_ts,
    forces_ts,
    last_images_positions,
):
    reactant, product, calc_defaults = setup_test_environment

    # for i in [reactant, product]:
    #    i.calc = NewtonNet(**calc_defaults)
    #    #i.calc = EMT()
    opt_defaults = {"optimizer": Sella, "optimizer_kwargs": ({"order": 0})}
    # opt_defaults = {"optimizer": BFGS}
    relax_job_kwargs = {"calc_kwargs": calc_defaults, "opt_kwargs": opt_defaults}
    optimized_r = strip_decorator(relax_job)(reactant, **relax_job_kwargs)["atoms"]
    optimized_p = strip_decorator(relax_job)(product, **relax_job_kwargs)["atoms"]
    # optimized_r = summarize_opt_run(run_opt(reactant, **opt_defaults))["atoms"]
    # optimized_p = summarize_opt_run(run_opt(product, **opt_defaults))["atoms"]
    # optimized_r.calc = EMT()
    # optimized_p.calc = EMT()
    optimized_r.calc = NewtonNet(**calc_defaults)
    optimized_p.calc = NewtonNet(**calc_defaults)

    images = _geodesic_interpolate_wrapper(
        optimized_r.copy(), optimized_p.copy(), nimages=n_intermediate
    )
    for image in images:
        image.calc = NewtonNet(**calc_defaults)
        # image.calc = EMT()
    assert 1 == 1
    # assert optimized_p.positions[0][1] == pytest.approx(
    # assert images[0].positions[0][1] == pytest.approx(first_image_positions, abs=1e-2)

    # assert optimized_p.get_potential_energy() == pytest.approx(first_image_pot_energy, abs=1e-2), "reactant pot. energy"

    # assert optimized_p.get_forces()[0, 1] == pytest.approx(first_image_forces, abs=1e-3), "reactant forces"

    # assert images[0].positions[0][1] == pytest.approx(
    #     last_images_positions,
    #     abs=1e-2,
    # )

    neb_kwargs = {"method": "aseneb", "precon": None}

    dyn = run_path_opt(images, optimizer=NEBOptimizer, neb_kwargs=neb_kwargs)

    neb_summary = summarize_path_opt_run(dyn)

    # assert neb_summary["trajectory_results"][1]["energy"] == pytest.approx(
    #    1.09889737,
    #    abs=0.1,
    # )
    assert neb_summary["trajectory_results"][1]["energy"] == pytest.approx(
        -24.650358983, abs=1
    )
    # assert images[1].positions[0][1] == pytest.approx(second_images_positions, abs=1e-1)

    # assert np.argmax(
    #    [image.get_potential_energy() for image in images]
    # ) == pytest.approx(index_ts), "Index of the transition state"

    # assert np.max([image.get_potential_energy() for image in images]) == pytest.approx(
    #    pot_energy_ts, abs=1), "Potential energy of the transition state"

    # assert images[
    #    np.argmax([image.get_potential_energy() for image in images])
    # ].get_forces()[0, 1] == pytest.approx(
    #    forces_ts, abs=1), "Force component in the transition state"

    # # Ensure the log file is correctly handled
    # if expected_logfile is None:
    #     assert logdir is None
    # else:
    #     assert logdir is not None
    #     log_filename = f"neb_band_{method}_{optimizer.__name__}_{precon}.txt"
    #     logfile_path = Path(logdir) / log_filename
    #     assert logfile_path.exists()
    #     # 'Could not find the optimization output file for NEB'
    #
    #     assert os.path.exists(
    #         f"{logdir}/optimized_path_{method}_{optimizer.__name__}_{precon}.xyz"
    #     ), "Could not find the xyz file for converged NEB calculation."
    # assert neb_summary[0] == 1


def test_run_calc(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    prep_files()

    with change_settings({"RESULTS_DIR": tmp_path}):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.1
        atoms.calc = EMT()

        new_atoms = run_calc(atoms, copy_files={Path(): "test_file.txt"})
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

        new_atoms = run_calc(atoms, copy_files={Path(): "test_file.txt"})
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

        dyn = run_opt(atoms, copy_files={Path(): "test_file.txt"})
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

    dyn = run_opt(
        atoms,
        optimizer=BFGS,
        copy_files={Path(): "test_file.txt"},
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None

    dyn = run_opt(
        traj[-1],
        optimizer=BFGSLineSearch,
        copy_files={Path(): "test_file.txt"},
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
    vib = run_vib(o2, copy_files={Path(): "test_file.txt"})
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
        run_calc(atoms, copy_files={Path(): "test_file.txt"})
    assert "Cannot find file" in caplog.text

    # No file again
    with caplog.at_level(logging.WARNING):
        run_opt(atoms, copy_files={Path(): "test_file.txt"})
    assert "Cannot find file" in caplog.text

    # No trajectory kwarg
    with pytest.raises(ValueError):
        run_opt(
            atoms,
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

        run_calc(atoms, copy_files={Path(): "test_file.txt"})
        results_dir = _find_results_dir()
        assert atoms.calc.results is not None
        assert not os.path.exists(os.path.join(results_dir, "test_file.txt"))
        assert os.path.exists(os.path.join(results_dir, "test_file.txt.gz"))

    with change_settings({"CREATE_UNIQUE_DIR": False, "RESULTS_DIR": tmp_path}):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.1
        atoms.calc = EMT()

        run_calc(atoms, copy_files={Path(): "test_file.txt"})
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
