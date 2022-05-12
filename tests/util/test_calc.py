import os
from shutil import rmtree

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones

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
    rmtree("blank_dir")


def test_run_calc():

    atoms = bulk("Cu")
    calc = EMT()
    atoms.calc = calc

    atoms = run_calc(
        atoms, scratch_dir="test_calc", gzip=False, copy_files=["test_file.txt"]
    )
    assert atoms.calc.results is not None
    assert os.path.exists("test_file.txt")
    assert not os.path.exists("test_file.txt.gz")

    atoms = run_calc(
        atoms, scratch_dir="test_calc", gzip=False, copy_files=["test_file.txt"]
    )
    assert atoms.calc.results is not None
    assert os.path.exists("test_file.txt")
    assert not os.path.exists("test_file.txt.gz")

    atoms = run_ase_opt(atoms, scratch_dir="test_calc", copy_files=["test_file.txt"])
    assert atoms[-1].calc.results is not None
    assert os.path.exists("test_file.txt")
    assert os.path.exists("test_file.txt.gz")
    os.remove("test_file.txt.gz")

    h2 = molecule("H2")
    h2.calc = LennardJones()
    vib = run_ase_vib(h2, scratch_dir="test_calc", copy_files=["test_file.txt"])
    assert np.real(vib.get_frequencies()[-1]) == pytest.approx(152229.50461546227)
    assert os.path.exists("test_file.txt")
    assert os.path.exists("test_file.txt.gz")
    os.remove("test_file.txt.gz")


def test_bad_run_calc():
    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        atoms = run_calc(atoms)
