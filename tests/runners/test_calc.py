import pytest

from quacc import SETTINGS

DEFAULT_SETTINGS = SETTINGS.copy()


def prep_files():
    import os

    # Make some test files to play with
    if not os.path.exists("test_calc"):
        os.mkdir("test_calc")
    with open("test_file.txt", "w") as f:
        f.write("test")


def teardown_function():
    import os
    from shutil import rmtree

    from quacc import SETTINGS

    if os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, "test_calc")):
        rmtree(os.path.join(SETTINGS.RESULTS_DIR, "test_calc"), ignore_errors=True)
    for f in ["test_file.txt", "test_file.txt.gz"]:
        if os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, f)):
            os.remove(os.path.join(SETTINGS.RESULTS_DIR, f))


def test_run_calc(tmpdir):
    import os

    import numpy as np
    from ase.build import bulk
    from ase.calculators.emt import EMT

    from quacc import SETTINGS
    from quacc.runners.calc import run_calc

    tmpdir.chdir()
    prep_files()

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    new_atoms = run_calc(atoms, copy_files=["test_file.txt"])
    assert atoms.calc.results is not None
    assert new_atoms.calc.results is not None
    assert not os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, "test_file.txt"))
    assert os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, "test_file.txt.gz"))
    assert np.array_equal(new_atoms.get_positions(), atoms.get_positions()) is True
    assert np.array_equal(new_atoms.cell.array, atoms.cell.array) is True


def test_run_calc_no_gzip(tmpdir):
    import os

    import numpy as np
    from ase.build import bulk
    from ase.calculators.emt import EMT

    from quacc import SETTINGS
    from quacc.runners.calc import run_calc

    tmpdir.chdir()
    prep_files()

    SETTINGS.GZIP_FILES = False

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    new_atoms = run_calc(atoms, copy_files=["test_file.txt"])
    assert atoms.calc.results is not None
    assert new_atoms.calc.results is not None
    assert os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, "test_file.txt"))
    assert not os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, "test_file.txt.gz"))
    assert np.array_equal(new_atoms.get_positions(), atoms.get_positions()) is True
    assert np.array_equal(new_atoms.cell.array, atoms.cell.array) is True
    SETTINGS.GZIP_FILES = DEFAULT_SETTINGS.GZIP_FILES


def test_run_ase_opt1(tmpdir):
    import os

    import numpy as np
    from ase.build import bulk
    from ase.calculators.emt import EMT

    from quacc import SETTINGS
    from quacc.runners.calc import run_ase_opt

    tmpdir.chdir()
    prep_files()

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    dyn = run_ase_opt(atoms, copy_files=["test_file.txt"])
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert not os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, "test_file.txt"))
    assert os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, "test_file.txt.gz"))
    assert np.array_equal(traj[-1].get_positions(), atoms.get_positions()) is False
    assert np.array_equal(traj[-1].cell.array, atoms.cell.array) is True


def test_run_ase_opt2(tmpdir):
    from ase.build import bulk
    from ase.calculators.emt import EMT
    from ase.optimize import BFGS, BFGSLineSearch

    from quacc.runners.calc import run_ase_opt

    tmpdir.chdir()
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    dyn = run_ase_opt(
        atoms,
        optimizer=BFGS,
        copy_files=["test_file.txt"],
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None

    dyn = run_ase_opt(
        traj[-1],
        optimizer=BFGSLineSearch,
        copy_files=["test_file.txt"],
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None


def test_run_ase_vib(tmpdir):
    import os

    import numpy as np
    from ase.build import molecule
    from ase.calculators.lj import LennardJones

    from quacc import SETTINGS
    from quacc.runners.calc import run_ase_vib

    tmpdir.chdir()
    prep_files()

    o2 = molecule("O2")
    o2.calc = LennardJones()
    vib = run_ase_vib(o2, copy_files=["test_file.txt"])
    assert np.real(vib.get_frequencies()[-1]) == pytest.approx(255.6863883406967)
    assert np.array_equal(vib.atoms.get_positions(), o2.get_positions()) is True
    assert not os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, "test_file.txt"))
    assert os.path.exists(os.path.join(SETTINGS.RESULTS_DIR, "test_file.txt.gz"))


def test_bad_runs(tmpdir):
    from ase.build import bulk
    from ase.calculators.emt import EMT
    from ase.optimize import BFGSLineSearch

    from quacc.runners.calc import run_ase_opt, run_calc

    tmpdir.chdir()

    atoms = bulk("Cu")

    # No calculator
    with pytest.raises(ValueError):
        run_calc(atoms)

    atoms.calc = EMT()

    # No file
    with pytest.warns(UserWarning):
        run_calc(atoms, copy_files=["test_file.txt"])

    # No file again
    with pytest.warns(UserWarning):
        run_ase_opt(atoms, copy_files=["test_file.txt"])

    # No trajectory kwarg
    with pytest.raises(ValueError):
        run_ase_opt(
            atoms,
            optimizer=BFGSLineSearch,
            optimizer_kwargs={
                "restart": None,
                "trajectory": "test.traj",
            },
        )


def test_unique_workdir(tmpdir):
    import numpy as np
    from ase.build import bulk, molecule
    from ase.calculators.emt import EMT
    from ase.calculators.lj import LennardJones

    from quacc import SETTINGS
    from quacc.runners.calc import run_ase_opt, run_ase_vib, run_calc

    SETTINGS.CREATE_UNIQUE_WORKDIR = True
    tmpdir.chdir()
    prep_files()

    # Static
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    run_calc(atoms, copy_files=["test_file.txt"])
    assert atoms.calc.results is not None

    # Opt
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    dyn = run_ase_opt(atoms, copy_files=["test_file.txt"])
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert np.array_equal(traj[-1].get_positions(), atoms.get_positions()) is False
    assert np.array_equal(traj[-1].cell.array, atoms.cell.array) is True

    # Vib
    o2 = molecule("O2")
    o2.calc = LennardJones()
    vib = run_ase_vib(o2)
    assert np.real(vib.get_frequencies()[-1]) == pytest.approx(255.6863883406967)
    assert np.array_equal(vib.atoms.get_positions(), o2.get_positions()) is True

    SETTINGS.CREATE_UNIQUE_WORKDIR = DEFAULT_SETTINGS.CREATE_UNIQUE_WORKDIR
