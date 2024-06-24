from __future__ import annotations

import os
from copy import deepcopy
from pathlib import Path

import pytest
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.io import read
from ase.optimize import BFGS
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.units import invcm
from ase.vibrations import Vibrations
from maggma.stores import MemoryStore
from monty.json import MontyDecoder, jsanitize
from monty.serialization import loadfn
from numpy.testing import assert_allclose as assert_close

from quacc.schemas.ase import (
    _summarize_harmonic_thermo,
    _summarize_ideal_gas_thermo,
    _summarize_vib_run,
    summarize_opt_run,
    summarize_run,
)

FILE_DIR = Path(__file__).parent

RUN1 = FILE_DIR / "test_files" / "vasp_run1"


def test_summarize_run(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    # Make sure metadata is made
    initial_atoms = read(os.path.join(RUN1, "POSCAR.gz"))
    atoms = read(os.path.join(RUN1, "OUTCAR.gz"))
    results = summarize_run(atoms, initial_atoms)
    assert results["nsites"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["results"]["energy"] == atoms.get_potential_energy()
    assert "pymatgen_version" in results["builder_meta"]
    assert results["input_atoms"]["atoms"] == initial_atoms
    assert Path(results["dir_name"]).is_dir()

    json_results = loadfn(Path(results["dir_name"], "quacc_results.json.gz"))
    assert json_results.keys() == results.keys()

    assert json_results["nsites"] == results["nsites"]
    assert json_results["results"]["energy"] == results["results"]["energy"]
    assert json_results["atoms"].info == results["atoms"].info


def test_summarize_run2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Test DB
    initial_atoms = read(os.path.join(RUN1, "POSCAR.gz"))
    atoms = read(os.path.join(RUN1, "OUTCAR.gz"))
    store = MemoryStore()
    summarize_run(atoms, initial_atoms, store=store)
    assert store.count() == 1


def test_summarize_run3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make sure info tags are handled appropriately
    initial_atoms = read(os.path.join(RUN1, "POSCAR.gz"))
    atoms = read(os.path.join(RUN1, "OUTCAR.gz"))
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = summarize_run(atoms, initial_atoms)
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}


def test_summarize_run4(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Make sure magnetic moments are handled appropriately
    initial_atoms = read(os.path.join(RUN1, "POSCAR.gz"))
    atoms = read(os.path.join(RUN1, "OUTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.calc.results["magmoms"] = [2.0] * len(atoms)
    results = summarize_run(atoms, initial_atoms, move_magmoms=True)

    assert atoms.calc is not None
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    assert results["atoms"].get_initial_magnetic_moments().tolist() == [2.0] * len(
        atoms
    )
    assert results["atoms"].calc is None


def test_summarize_run5(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Make sure Atoms magmoms were not moved if specified
    initial_atoms = read(os.path.join(RUN1, "POSCAR.gz"))
    atoms = read(os.path.join(RUN1, "OUTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    results = summarize_run(
        atoms, initial_atoms, move_magmoms=False, additional_fields={"test": "hi"}
    )

    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    assert results["atoms"].get_initial_magnetic_moments().tolist() == [3.14] * len(
        atoms
    )
    assert results.get("test") == "hi"

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_summarize_opt_run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make sure metadata is made
    atoms = bulk("Cu") * (2, 2, 1)
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    dyn = BFGS(atoms, trajectory="test.traj")
    dyn.run(steps=100)
    traj = read("test.traj", index=":")

    results = summarize_opt_run(dyn)
    assert results["nsites"] == len(atoms)
    assert results["atoms"] == traj[-1]
    assert results["results"]["energy"] == atoms.get_potential_energy()
    assert len(results["trajectory"]) == len(traj)
    assert len(results["trajectory_results"]) == len(traj)
    assert results["trajectory_results"][-1]["energy"] == results["results"]["energy"]
    assert "nid" in results
    assert "dir_name" in results
    assert "pymatgen_version" in results["builder_meta"]
    assert results["parameters_opt"]["fmax"] == dyn.fmax
    assert results["parameters_opt"]["max_steps"] == 100

    json_results = loadfn(Path(results["dir_name"], "quacc_results.json.gz"))

    assert json_results.keys() == results.keys()

    # assert things on the trajectory are the same
    assert json_results["trajectory"] == results["trajectory"]
    assert (
        json_results["trajectory_results"][-1]["energy"]
        == results["trajectory_results"][-1]["energy"]
    )

    # Test DB
    atoms = bulk("Cu") * (2, 2, 1)
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    dyn = BFGS(atoms, trajectory="test.traj")
    dyn.run(steps=5)
    traj = read("test.traj", index=":")

    store = MemoryStore()
    summarize_opt_run(dyn, store=store, check_convergence=False)
    assert store.count() == 1

    # Test no convergence
    atoms = bulk("Cu") * (2, 2, 1)
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    dyn = BFGS(atoms, trajectory="test.traj")
    dyn.run(steps=5)
    traj = read("test.traj", index=":")

    with pytest.raises(RuntimeError, match="Optimization did not converge"):
        summarize_opt_run(dyn)

    # Make sure info tags are handled appropriately
    atoms = bulk("Cu") * (2, 2, 1)
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    dyn = BFGS(atoms, trajectory="test.traj")
    dyn.run()

    results = summarize_opt_run(dyn)
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)

    with pytest.raises(FileNotFoundError):
        dyn.trajectory.filename = "not_a_file.traj"
        summarize_opt_run(dyn)


def test_summarize_vib_run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make sure metadata is made
    atoms = molecule("N2")
    atoms.calc = EMT()
    input_atoms = deepcopy(atoms)
    vib = Vibrations(atoms)
    vib.run()

    results = _summarize_vib_run(vib)
    assert results["atoms"] == input_atoms
    assert results["natoms"] == len(atoms)
    assert results["parameters_vib"]["delta"] == vib.delta
    assert results["parameters_vib"]["direction"] == "central"
    assert results["parameters_vib"]["method"] == "standard"
    assert results["parameters_vib"]["ndof"] == 6
    assert results["parameters_vib"]["nfree"] == 2
    assert "nid" in results
    assert "dir_name" in results
    assert "pymatgen_version" in results["builder_meta"]
    assert len(results["results"]["vib_freqs_raw"]) == 6
    assert results["results"]["vib_freqs_raw"][0] == pytest.approx(0, rel=1e-5)
    assert results["results"]["vib_freqs_raw"][-1] == pytest.approx(928.1447554058556)
    assert len(results["results"]["vib_energies_raw"]) == 6
    assert results["results"]["vib_energies_raw"][0] == pytest.approx(0, rel=1e-5)
    assert results["results"]["vib_energies_raw"][-1] == pytest.approx(
        0.11507528256667966
    )
    assert results["results"]["n_imag"] == 0
    assert results["results"]["imag_vib_freqs"] == []
    assert len(results["results"]["vib_freqs"]) == 1
    assert results["results"]["vib_freqs"][0] == pytest.approx(928.1447554058556)
    assert len(results["results"]["vib_energies"]) == 1
    assert results["results"]["vib_energies"][0] == pytest.approx(0.11507528256667966)

    # Test DB
    atoms = molecule("N2")
    atoms.calc = EMT()
    input_atoms = deepcopy(atoms)
    vib = Vibrations(atoms)
    vib.run()

    # Make sure info tags are handled appropriately
    atoms = molecule("N2")
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms.calc = EMT()
    vib = Vibrations(atoms)
    vib.run()

    results = _summarize_vib_run(vib)
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

    results = _summarize_vib_run(vib)
    assert results["atoms"] == input_atoms
    assert results["nsites"] == len(atoms)
    assert results["parameters_vib"]["delta"] == vib.delta
    assert len(results["results"]["vib_freqs_raw"]) == 6
    assert len(results["results"]["vib_energies_raw"]) == 6
    assert len(results["results"]["vib_freqs"]) == 6
    assert len(results["results"]["vib_energies"]) == 6


def test_summarize_ideal_gas_thermo(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make sure metadata is made
    atoms = molecule("N2")
    igt = IdealGasThermo([0.34], "linear", atoms=atoms, spin=0, symmetrynumber=2)
    results = _summarize_ideal_gas_thermo(igt)
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["parameters_thermo"]["vib_energies"] == [0.34]
    assert results["parameters_thermo"]["vib_freqs"] == [0.34 / invcm]
    assert results["results"]["energy"] == 0
    assert "pymatgen_version" in results["builder_meta"]

    # Make sure right number of vib energies are reported
    atoms = molecule("N2")
    igt = IdealGasThermo(
        [0.0, 0.34], "linear", atoms=atoms, potentialenergy=-1, spin=0, symmetrynumber=2
    )
    results = _summarize_ideal_gas_thermo(igt)
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["parameters_thermo"]["vib_energies"] == [0.34]
    assert results["parameters_thermo"]["vib_freqs"] == [0.34 / invcm]
    assert results["results"]["energy"] == -1

    # # Make sure info tags are handled appropriately
    atoms = molecule("N2")
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms.calc = EMT()
    igt = IdealGasThermo(
        [0.0, 0.34], "linear", atoms=atoms, potentialenergy=-1, spin=0, symmetrynumber=2
    )
    results = _summarize_ideal_gas_thermo(igt)
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
    results = _summarize_ideal_gas_thermo(igt, temperature=1000.0, pressure=20.0)
    assert results["natoms"] == len(atoms)
    assert results["atoms"] == atoms
    assert len(results["parameters_thermo"]["vib_energies"]) == 6
    assert results["parameters_thermo"]["vib_energies"][0] == vib_energies[-6]
    assert results["parameters_thermo"]["vib_energies"][-1] == vib_energies[-1]
    assert results["results"]["energy"] == -10.0
    assert results["results"]["enthalpy"] == pytest.approx(-8.749341973959462)
    assert results["results"]["entropy"] == pytest.approx(0.0023506788982171896)
    assert results["results"]["gibbs_energy"] == pytest.approx(-11.100020872176652)
    assert results["parameters_thermo"]["temperature"] == 1000.0
    assert results["parameters_thermo"]["pressure"] == 20.0
    assert results["parameters_thermo"]["sigma"] == 6
    assert results["parameters_thermo"]["spin_multiplicity"] == 2

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)

    with pytest.raises(ValueError):
        _summarize_ideal_gas_thermo(igt, charge_and_multiplicity=[0, 1])


def test_summarize_harmonic_thermo(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make sure metadata is made
    ht = HarmonicThermo([0.34])
    results = _summarize_harmonic_thermo(ht)
    assert results["parameters_thermo"]["vib_energies"] == [0.34]
    assert results["parameters_thermo"]["vib_freqs"] == [0.34 / invcm]
    assert results["results"]["energy"] == 0
    assert_close(results["results"]["helmholtz_energy"], 0.16999995401497991, rtol=1e-5)
    assert_close(results["results"]["internal_energy"], 0.1700006085385999, rtol=1e-5)
    assert_close(results["results"]["entropy"], 2.1952829783392438e-09, rtol=1e-5)
    assert_close(results["results"]["zpe"], 0.17, rtol=1e-5)

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_errors(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        summarize_run(atoms, atoms)

    initial_atoms = read(os.path.join(RUN1, "POSCAR.gz"))
    atoms = read(os.path.join(RUN1, "OUTCAR.gz"))
    atoms.calc.results = {}
    with pytest.raises(ValueError):
        summarize_run(atoms, initial_atoms)
