import os
from shutil import rmtree

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones

from quacc.util.calc import ideal_gas_thermo, run_ase_opt, run_ase_vib, run_calc

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
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    calc = EMT()
    atoms.calc = calc

    new_atoms = run_calc(
        atoms, scratch_dir="test_calc", gzip=False, copy_files=["test_file.txt"]
    )
    assert atoms.calc.results is not None
    assert os.path.exists("test_file.txt")
    assert not os.path.exists("test_file.txt.gz")
    assert np.array_equal(new_atoms.get_positions(), atoms.get_positions()) is True
    assert np.array_equal(new_atoms.cell.array, atoms.cell.array) is True

    new_atoms = run_calc(
        atoms, scratch_dir="test_calc", gzip=False, copy_files=["test_file.txt"]
    )
    assert new_atoms.calc.results is not None
    assert os.path.exists("test_file.txt")
    assert not os.path.exists("test_file.txt.gz")
    assert np.array_equal(new_atoms.get_positions(), atoms.get_positions()) is True
    assert np.array_equal(new_atoms.cell.array, atoms.cell.array) is True

    new_atoms = run_ase_opt(
        atoms, scratch_dir="test_calc", copy_files=["test_file.txt"]
    )
    assert new_atoms[-1].calc.results is not None
    assert os.path.exists("test_file.txt")
    assert os.path.exists("test_file.txt.gz")
    assert np.array_equal(new_atoms[-1].get_positions(), atoms.get_positions()) is False
    assert np.array_equal(new_atoms[-1].cell.array, atoms.cell.array) is True
    os.remove("test_file.txt.gz")

    o2 = molecule("O2")
    o2.calc = LennardJones()
    vib = run_ase_vib(o2, scratch_dir="test_calc", copy_files=["test_file.txt"])
    assert np.real(vib.get_frequencies()[-1]) == pytest.approx(255.6863883406967)
    assert np.array_equal(vib.atoms.get_positions(), o2.get_positions()) is True
    assert os.path.exists("test_file.txt")
    assert os.path.exists("test_file.txt.gz")
    os.remove("test_file.txt.gz")


def test_bad_run_calc():
    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        atoms = run_calc(atoms)


# def test_ideal_gas_thermo():
#     # Note: More detailed tests are in the test for the xtb ThermoJob
#     atoms = molecule("CH4")
#     dummy_freqs = [
#         0,
#         0,
#         0,
#         0,
#         0,
#         10j,
#         200j,
#         500 + 0j,
#         1000 + 0j,
#         1500,
#         2000,
#         2500,
#         3000,
#         3500,
#         4000,
#     ]
#     igt = ideal_gas_thermo(atoms, dummy_freqs)
#     assert igt["results"]["frequencies"] == [
#         0,
#         0,
#         0,
#         0,
#         0,
#         -10,
#         -200,
#         500,
#         1000,
#         1500,
#         2000,
#         2500,
#         3000,
#         3500,
#         4000,
#     ]
#     assert igt["results"]["true_frequencies"] == [
#         -200,
#         500,
#         1000,
#         1500,
#         2000,
#         2500,
#         3000,
#         3500,
#         4000,
#     ]
#     assert igt["results"]["n_imag"] == 1  # 200j
#     assert igt["symmetry"]["linear"] is False
#     assert igt["symmetry"]["point_group"] == "Td"
#     assert igt["symmetry"]["rotation_number"] == 12

#     atoms = molecule("CH4")
#     dummy_freqs = [
#         0,
#         0,
#         0,
#         0,
#         0,
#         -10,
#         0.0 + 200j,
#         500 + 0j,
#         1000 + 0j,
#         1500,
#         2000,
#         2500,
#         3000,
#         3500,
#         4000,
#     ]
#     igt = ideal_gas_thermo(atoms, dummy_freqs)
#     assert igt["results"]["frequencies"] == [
#         0,
#         0,
#         0,
#         0,
#         0,
#         -10,
#         -200,
#         500,
#         1000,
#         1500,
#         2000,
#         2500,
#         3000,
#         3500,
#         4000,
#     ]
#     assert igt["results"]["true_frequencies"] == [
#         -200,
#         500,
#         1000,
#         1500,
#         2000,
#         2500,
#         3000,
#         3500,
#         4000,
#     ]
#     assert igt["results"]["n_imag"] == 1  # 200j
