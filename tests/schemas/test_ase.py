import os
from pathlib import Path
from shutil import rmtree

import pytest
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.io import read
from ase.optimize import BFGS
from ase.thermochemistry import IdealGasThermo
from ase.units import invcm
from ase.vibrations import Vibrations
from monty.json import MontyDecoder, jsanitize

from quacc.schemas.ase import (
    summarize_opt_run,
    summarize_run,
    summarize_thermo_run,
    summarize_vib_run,
)

FILE_DIR = Path(__file__).resolve().parent

run1 = os.path.join(FILE_DIR, "vasp_run1")


def teardown_module():
    for f in os.listdir("."):
        if ".log" in f or ".pckl" in f or ".traj" in f or ".gz" in f:
            os.remove(f)
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir" or f == "vib":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_summarize_run():
    # Make sure metadata is made
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms)
    assert results["nsites"] == len(atoms)
    assert results["atoms"] == atoms

    # Make sure initial atoms object is stored if specified
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms, atoms)
    assert results["nsites"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["input_structure"]["atoms"] == atoms

    # Make sure info tags are handled appropriately
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = summarize_run(atoms)
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # Make sure magnetic moments are handled appropriately
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.calc.results["magmoms"] = [2.0] * len(atoms)
    results = summarize_run(atoms)

    assert atoms.calc is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    assert results["atoms"].get_initial_magnetic_moments().tolist() == [2.0] * len(
        atoms
    )
    assert results["atoms"].calc is None

    # Make sure Atoms magmoms were not moved if specified
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = summarize_run(atoms, prep_next_run=False)
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    assert results["atoms"].get_initial_magnetic_moments().tolist() == [3.14] * len(
        atoms
    )

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_summarize_opt_run():
    # Make sure metadata is made
    atoms = bulk("Cu") * (2, 2, 1)
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    BFGS(atoms, trajectory="test.traj").run()
    traj = read("test.traj", index=":")

    results = summarize_opt_run(traj)
    assert results["nsites"] == len(atoms)
    assert results["atoms"] == traj[-1]
    assert len(results["trajectory"]) == len(traj)
    assert len(results["trajectory_results"]) == len(traj)
    assert results["trajectory_results"][-1] == results["results"]

    # Make sure info tags are handled appropriately
    atoms = bulk("Cu") * (2, 2, 1)
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    BFGS(atoms, trajectory="test.traj").run()
    traj = read("test.traj", index=":")

    results = summarize_opt_run(traj)
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # Make sure magnetic moments are handled appropriately
    atoms = bulk("Cu") * (2, 2, 1)
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    BFGS(atoms, trajectory="test.traj").run()
    traj = read("test.traj", index=":")
    traj[-1].calc.results["magmoms"] = [2.0] * len(atoms)

    results = summarize_opt_run(traj)

    assert atoms.calc is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    assert results["atoms"].get_initial_magnetic_moments().tolist() == [2.0] * len(
        atoms
    )
    assert results["atoms"].calc is None

    # Make sure Atoms magmoms were not moved if specified
    atoms = bulk("Cu") * (2, 2, 1)
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    BFGS(atoms, trajectory="test.traj").run()
    traj = read("test.traj", index=":")
    traj[-1].calc.results["magmoms"] = [2.0] * len(atoms)

    results = summarize_opt_run(traj, prep_next_run=False)
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    assert results["atoms"].get_initial_magnetic_moments().tolist() == [3.14] * len(
        atoms
    )

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_summarize_vib_run():
    # Make sure metadata is made
    atoms = molecule("N2")
    atoms.calc = EMT()
    vib = Vibrations(atoms)
    vib.run()

    results = summarize_vib_run(vib)
    assert results["natoms"] == len(atoms)
    assert results["parameters"]["delta"] == vib.delta

    # Make sure info tags are handled appropriately
    atoms = molecule("N2")
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms.calc = EMT()
    vib = Vibrations(atoms)
    vib.run()

    results = summarize_vib_run(vib)
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_summarize_thermo_run():
    # Make sure metadata is made
    atoms = molecule("N2")
    igt = IdealGasThermo([0.34], "linear", atoms=atoms, spin=0, symmetrynumber=2)
    results = summarize_thermo_run(igt)
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["results"]["vib_energies"] == [0.34]
    assert results["results"]["vib_freqs"] == [0.34 / invcm]
    assert results["results"]["energy"] == 0

    # Make sure right number of vib energies are reported
    atoms = molecule("N2")
    igt = IdealGasThermo(
        [0.0, 0.34], "linear", atoms=atoms, potentialenergy=-1, spin=0, symmetrynumber=2
    )
    results = summarize_thermo_run(igt)
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["results"]["vib_energies"] == [0.34]
    assert results["results"]["vib_freqs"] == [0.34 / invcm]
    assert results["results"]["energy"] == -1

    # # Make sure info tags are handled appropriately
    atoms = molecule("N2")
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms.calc = EMT()
    igt = IdealGasThermo(
        [0.0, 0.34], "linear", atoms=atoms, potentialenergy=-1, spin=0, symmetrynumber=2
    )
    results = summarize_thermo_run(igt)
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_errors():
    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        summarize_run(atoms)

    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    atoms.calc.results = {}
    with pytest.raises(ValueError):
        summarize_run(atoms)
