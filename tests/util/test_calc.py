import os
from shutil import rmtree

import pytest
from ase.build import bulk

from quacc.calculators.vasp import SmartVasp
from quacc.util.calc import run_calc

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
    atoms = SmartVasp(atoms)

    atoms = run_calc(atoms, scratch_dir="test_calc", copy_from_store_dir=True)
    assert atoms.calc.results is not None
    assert os.path.exists("test_file.txt")
    assert os.path.exists("test_file.txt.gz")
    os.remove("test_file.txt.gz")

    atoms = run_calc(
        atoms, scratch_dir="test_calc", gzip=False, copy_from_store_dir=True
    )
    assert atoms.calc.results is not None
    assert os.path.exists("test_file.txt")
    assert not os.path.exists("test_file.txt.gz")


def test_bad_run_calc(monkeypatch):
    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        atoms = run_calc(atoms)
