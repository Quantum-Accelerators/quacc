import os
from shutil import rmtree

import pytest
from ase.build import bulk

from quacc.calculators.vasp import SmartVasp
from quacc.util.calc import run_calc


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
    os.chdir("..")
    rmtree("blank_dir")


def test_run_calc():

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms)

    atoms = run_calc(atoms)
    if atoms.calc.results is None:
        raise AssertionError
    if not os.path.exists("test_file.txt"):
        raise AssertionError
    if os.path.exists("test_file.txt.gz"):
        raise AssertionError

    atoms = run_calc(atoms, gzip=True)
    if atoms.calc.results is None:
        raise AssertionError
    if not os.path.exists("test_file.txt"):
        raise AssertionError
    if not os.path.exists("test_file.txt.gz"):
        raise AssertionError

    atoms = run_calc(atoms, store_dir=".", scratch_dir="test_calc", gzip=True)
    if atoms.calc.results is None:
        raise AssertionError
    if not os.path.exists("test_file.txt"):
        raise AssertionError
    if not os.path.exists("test_file.txt.gz"):
        raise AssertionError


def test_bad_run_calc(monkeypatch):
    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        atoms = run_calc(atoms)

    atoms = SmartVasp(atoms)

    monkeypatch.setenv("SCRATCH", "nonexistant_dir")
    with pytest.raises(OSError):
        atoms = run_calc(atoms)
