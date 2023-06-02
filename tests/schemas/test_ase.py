import os
from copy import deepcopy
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
    assert results["results"]["energy"] == atoms.get_potential_energy()
    assert "pull_request" in results["builder_meta"]

    # Test remove_empties
    atoms = read(os.path.join(run1, "OUTCAR.gz"))
    results = summarize_run(atoms, remove_empties=True)
    assert results["nsites"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["results"]["energy"] == atoms.get_potential_energy()
    assert "pull_request" not in results["builder_meta"]

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
    dyn = BFGS(atoms, trajectory="test.traj")
    dyn.trajectory.filename = "test.traj"  # can remove after ASE MR 2901
    dyn.run()
    traj = read("test.traj", index=":")

    results = summarize_opt_run(dyn, check_convergence=True)
    assert results["nsites"] == len(atoms)
    assert results["atoms"] == traj[-1]
    assert results["results"]["energy"] == atoms.get_potential_energy()
    assert len(results["trajectory"]) == len(traj)
    assert len(results["trajectory_results"]) == len(traj)
    assert results["trajectory_results"][-1]["energy"] == results["results"]["energy"]
    assert "nid" in results
    assert "dir_name" in results
    assert "pull_request" in results["builder_meta"]

    # Test remove_empties
    atoms = bulk("Cu") * (2, 2, 1)
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    dyn = BFGS(atoms, trajectory="test.traj")
    dyn.trajectory.filename = "test.traj"  # can remove after ASE MR 2901
    dyn.run()
    traj = read("test.traj", index=":")

    results = summarize_opt_run(dyn, remove_empties=True)
    assert results["nsites"] == len(atoms)
    assert results["atoms"] == traj[-1]
    assert results["results"]["energy"] == atoms.get_potential_energy()
    assert len(results["trajectory"]) == len(traj)
    assert len(results["trajectory_results"]) == len(traj)
    assert results["trajectory_results"][-1]["energy"] == results["results"]["energy"]
    assert "nid" in results
    assert "dir_name" in results
    assert "pull_request" not in results["builder_meta"]

    # Make sure info tags are handled appropriately
    atoms = bulk("Cu") * (2, 2, 1)
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    dyn = BFGS(atoms, trajectory="test.traj")
    dyn.run()
    traj = read("test.traj", index=":")

    results = summarize_opt_run(dyn)
    assert results.get("atoms_info", {}) != {}
    assert results["atoms_info"].get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)

    with pytest.raises(FileNotFoundError):
        dyn.trajectory.filename = "not_a_file.traj"
        summarize_opt_run(dyn)


def test_summarize_vib_run():
    # Make sure metadata is made
    atoms = molecule("N2")
    atoms.calc = EMT()
    input_atoms = deepcopy(atoms)
    vib = Vibrations(atoms)
    vib.run()

    results = summarize_vib_run(vib)
    assert results["atoms"] == input_atoms
    assert results["natoms"] == len(atoms)
    assert results["parameters"]["delta"] == vib.delta
    assert results["parameters"]["direction"] == "central"
    assert results["parameters"]["method"] == "standard"
    assert results["parameters"]["ndof"] == 6
    assert results["parameters"]["nfree"] == 2
    assert "nid" in results
    assert "dir_name" in results
    assert "pull_request" in results["builder_meta"]
    assert len(results["results"]["vib_freqs_raw"]) == 6
    assert results["results"]["vib_freqs_raw"][0] == pytest.approx(
        -3.054403266390365e-06
    )
    assert results["results"]["vib_freqs_raw"][-1] == pytest.approx(928.1447554058556)
    assert len(results["results"]["vib_energies_raw"]) == 6
    assert results["results"]["vib_energies_raw"][0] == pytest.approx(
        -3.786977375083739e-10
    )
    assert results["results"]["vib_energies_raw"][-1] == pytest.approx(
        0.11507528256667966
    )
    assert results["results"]["n_imag"] == 0
    assert results["results"]["imag_vib_freqs"] == []
    assert len(results["results"]["vib_freqs"]) == 1
    assert results["results"]["vib_freqs"][0] == pytest.approx(928.1447554058556)
    assert len(results["results"]["vib_energies"]) == 1
    assert results["results"]["vib_energies"][0] == pytest.approx(0.11507528256667966)

    # Test remove_empties
    atoms = molecule("N2")
    atoms.calc = EMT()
    input_atoms = deepcopy(atoms)
    vib = Vibrations(atoms)
    vib.run()

    results = summarize_vib_run(vib, remove_empties=True)
    assert results["atoms"] == input_atoms
    assert results["natoms"] == len(atoms)
    assert results["parameters"]["delta"] == vib.delta
    assert results["parameters"]["direction"] == "central"
    assert results["parameters"]["method"] == "standard"
    assert results["parameters"]["ndof"] == 6
    assert results["parameters"]["nfree"] == 2
    assert "nid" in results
    assert "dir_name" in results
    assert "pull_request" not in results["builder_meta"]

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

    # Test a solid
    atoms = bulk("Cu") * (2, 1, 1)
    atoms.calc = EMT()
    input_atoms = deepcopy(atoms)
    vib = Vibrations(atoms)
    vib.run()

    results = summarize_vib_run(vib)
    assert results["atoms"] == input_atoms
    assert results["nsites"] == len(atoms)
    assert results["parameters"]["delta"] == vib.delta
    assert len(results["results"]["vib_freqs_raw"]) == 6
    assert len(results["results"]["vib_energies_raw"]) == 6
    assert len(results["results"]["vib_freqs"]) == 6
    assert len(results["results"]["vib_energies"]) == 6


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
    assert "pull_request" in results["builder_meta"]

    # Test remove_empties
    atoms = molecule("N2")
    igt = IdealGasThermo([0.34], "linear", atoms=atoms, spin=0, symmetrynumber=2)
    results = summarize_thermo_run(igt, remove_empties=True)
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["results"]["vib_energies"] == [0.34]
    assert results["results"]["vib_freqs"] == [0.34 / invcm]
    assert results["results"]["energy"] == 0
    assert "pull_request" not in results["builder_meta"]

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

    # Make sure spin works right
    atoms = molecule("CH3")
    vib_energies = [
        9.551077150221621e-06j,
        3.1825877476455407e-06j,
        2.7223332245579342e-06j,
        (0.03857802457526743 + 0j),
        (0.038762952842240087 + 0j),
        (0.03876411386029907 + 0j),
        (0.07135067701372912 + 0j),
        (0.1699785717790056 + 0j),
        (0.1700229358789492 + 0j),
        (0.3768719400148424 + 0j),
        (0.38803854931751625 + 0j),
        (0.3880868821616261 + 0j),
    ]
    igt = IdealGasThermo(
        vib_energies,
        "nonlinear",
        potentialenergy=-10.0,
        atoms=atoms,
        spin=0.5,
        symmetrynumber=6,
    )
    results = summarize_thermo_run(igt, temperature=1000.0, pressure=20.0)
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert len(results["results"]["vib_energies"]) == 6
    assert results["results"]["vib_energies"][0] == vib_energies[-6]
    assert results["results"]["vib_energies"][-1] == vib_energies[-1]
    assert results["results"]["energy"] == -10.0
    assert results["results"]["enthalpy"] == pytest.approx(-8.749341973959462)
    assert results["results"]["entropy"] == pytest.approx(0.0023506788982171896)
    assert results["results"]["gibbs_energy"] == pytest.approx(-11.100020872176652)
    assert results["parameters"]["temperature"] == 1000.0
    assert results["parameters"]["pressure"] == 20.0
    assert results["parameters"]["sigma"] == 6
    assert results["parameters"]["spin_multiplicity"] == 2
    assert "nid" in results
    assert "dir_name" in results

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
