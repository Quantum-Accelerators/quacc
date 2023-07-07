import os
from shutil import rmtree
from subprocess import Popen

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS, BFGSLineSearch
from monty.os.path import zpath

from quacc.util.calc import run_ase_opt, run_ase_vib, run_calc

CWD = os.getcwd()


def setup_module():
    # Run this test from a fresh directory
    if not os.path.exists("blank_dir"):
        os.mkdir("blank_dir")
    os.chdir("blank_dir")

    # Make some test files to play with
    if not os.path.exists("test_calc"):
        os.mkdir("test_calc")
    with open("test_file.txt", "w") as f:
        f.write("test")


def teardown_module():
    # Clean up
    os.chdir(CWD)
    for f in os.listdir("."):
        if ".log" in f or ".pckl" in f or ".traj" in f:
            os.remove(f)
    for f in os.listdir(CWD):
        if (
            "quacc-tmp" in f
            or "quacc_" in f
            or f == "tmp_dir"
            or f == "vib"
            or f == "blank_dir"
        ):
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_run_calc():
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    new_atoms = run_calc(
        atoms,
        scratch_dir="test_calc",
        gzip=False,
        copy_files=["test_file.txt"],
    )
    assert atoms.calc.results is not None
    assert os.path.exists("test_file.txt")
    assert not os.path.exists("test_file.txt.gz")
    assert np.array_equal(new_atoms.get_positions(), atoms.get_positions()) is True
    assert np.array_equal(new_atoms.cell.array, atoms.cell.array) is True

    new_atoms = run_calc(
        atoms,
        scratch_dir="test_calc",
        gzip=False,
        copy_files=["test_file.txt"],
    )
    assert new_atoms.calc.results is not None
    assert os.path.exists("test_file.txt")
    assert not os.path.exists("test_file.txt.gz")
    assert np.array_equal(new_atoms.get_positions(), atoms.get_positions()) is True
    assert np.array_equal(new_atoms.cell.array, atoms.cell.array) is True

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    new_atoms = run_calc(
        atoms,
        scratch_dir="new_test_calc",
        gzip=False,
        copy_files=["test_file.txt"],
    )
    assert atoms.calc.results is not None

    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        run_calc(
            atoms,
            scratch_dir="test_calc",
            gzip=False,
            copy_files=["test_file.txt"],
        )


def test_run_ase_opt():
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    dyn = run_ase_opt(atoms, scratch_dir="test_calc", copy_files=["test_file.txt"])
    Popen(f"gunzip {dyn.trajectory.filename}.gz", shell=True).wait()
    traj = read(dyn.trajectory.filename, index=":")
    assert traj[-1].calc.results is not None
    assert os.path.exists("test_file.txt")
    assert os.path.exists("test_file.txt.gz")
    assert np.array_equal(traj[-1].get_positions(), atoms.get_positions()) is False
    assert np.array_equal(traj[-1].cell.array, atoms.cell.array) is True
    os.remove("test_file.txt.gz")

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()

    dyn = run_ase_opt(
        atoms,
        optimizer=BFGS,
        scratch_dir="new_test_calc2",
        gzip=False,
        copy_files=["test_file.txt"],
        optimizer_kwargs={"restart": None},
    )
    Popen(f"gunzip {dyn.trajectory.filename}.gz", shell=True).wait()
    traj = read(dyn.trajectory.filename, index=":")
    assert traj[-1].calc.results is not None

    dyn = run_ase_opt(
        traj[-1],
        optimizer=BFGSLineSearch,
        scratch_dir="test_calc",
        gzip=False,
        copy_files=["test_file.txt"],
        optimizer_kwargs={"restart": None, "trajectory": "new_test.traj"},
    )
    Popen(f"gunzip {dyn.trajectory.filename}.gz", shell=True).wait()
    traj = read(dyn.trajectory.filename, index=":")
    assert traj[-1].calc.results is not None

    with pytest.raises(ValueError):
        dyn = run_ase_opt(
            traj[-1],
            optimizer=BFGSLineSearch,
            scratch_dir="test_calc",
            gzip=False,
            copy_files=["test_file.txt"],
            optimizer_kwargs={
                "restart": None,
                "trajectory": Trajectory("new_test2.traj", "w", atoms=traj[-1]),
            },
        )

    with pytest.raises(ValueError):
        run_ase_opt(
            bulk("Cu"),
            scratch_dir="test_calc",
            copy_files=["test_file.txt"],
        )


def test_run_ase_vib():
    o2 = molecule("O2")
    o2.calc = LennardJones()
    vib = run_ase_vib(o2, scratch_dir="test_calc_vib", copy_files=["test_file.txt"])
    assert np.real(vib.get_frequencies()[-1]) == pytest.approx(255.6863883406967)
    assert np.array_equal(vib.atoms.get_positions(), o2.get_positions()) is True
    assert os.path.exists("test_file.txt")
    assert os.path.exists("test_file.txt.gz")
    os.remove("test_file.txt.gz")

    with pytest.raises(ValueError):
        run_ase_vib(bulk("Cu"))


def test_bad_run_calc():
    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        atoms = run_calc(atoms)
