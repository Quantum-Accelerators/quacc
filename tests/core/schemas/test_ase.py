from __future__ import annotations

import os
from copy import deepcopy
from pathlib import Path

import pytest
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.io import read
from ase.mep import NEB
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from monty.json import MontyDecoder, jsanitize

from quacc.schemas.ase import Summarize, VibSummarize

FILE_DIR = Path(__file__).parent

RUN1 = FILE_DIR / "test_files" / "vasp_run1"


def test_summarize_run(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    # Make sure metadata is made
    initial_atoms = read(os.path.join(RUN1, "POSCAR.gz"))
    atoms = read(os.path.join(RUN1, "OUTCAR.gz"))
    results = Summarize().run(atoms, initial_atoms)
    assert results["structure_metadata"]["nsites"] == len(atoms)
    assert results["atoms"] == atoms
    assert results["results"]["energy"] == atoms.get_potential_energy()
    assert results["input_atoms"]["atoms"] == initial_atoms
    assert Path(results["dir_name"]).is_dir()


def test_summarize_run3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make sure info tags are handled appropriately
    initial_atoms = read(os.path.join(RUN1, "POSCAR.gz"))
    atoms = read(os.path.join(RUN1, "OUTCAR.gz"))
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    results = Summarize().run(atoms, initial_atoms)
    assert atoms.info.get("test_dict", None) == {"hi": "there", "foo": "bar"}
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}


def test_summarize_run4(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Make sure magnetic moments are handled appropriately
    initial_atoms = read(os.path.join(RUN1, "POSCAR.gz"))
    atoms = read(os.path.join(RUN1, "OUTCAR.gz"))
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms.calc.results["magmoms"] = [2.0] * len(atoms)
    results = Summarize(move_magmoms=True).run(atoms, initial_atoms)

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
    results = Summarize(move_magmoms=False, additional_fields={"test": "hi"}).run(
        atoms, initial_atoms
    )

    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)
    assert results["atoms"].get_initial_magnetic_moments().tolist() == [3.14] * len(
        atoms
    )
    assert results.get("test") == "hi"

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_summarize_opt_run1(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make sure metadata is made
    atoms = bulk("Cu") * (2, 2, 1)
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    dyn = BFGS(atoms, trajectory="test.traj")
    dyn.run(steps=100)
    traj = read("test.traj", index=":")

    results = Summarize().opt(dyn)
    assert results["structure_metadata"]["nsites"] == len(atoms)
    assert results["atoms"] == traj[-1]
    assert results["results"]["energy"] == atoms.get_potential_energy()
    assert len(results["trajectory"]) == len(traj)
    assert len(results["trajectory_results"]) == len(traj)
    assert results["trajectory_results"][-1]["energy"] == results["results"]["energy"]
    assert "nid" in results
    assert "dir_name" in results
    assert results["parameters_opt"]["fmax"] == dyn.fmax
    assert results["parameters_opt"]["max_steps"] == 100

    # Test custom traj
    assert (
        Summarize().opt(dyn, trajectory=traj, check_convergence=False)["trajectory"]
        == traj
    )


def test_summarize_opt_run3(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Test no convergence
    atoms = bulk("Cu") * (2, 2, 1)
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    dyn = BFGS(atoms, trajectory="test.traj")
    dyn.run(steps=5)

    with pytest.raises(RuntimeError, match="Optimization did not converge"):
        Summarize().opt(dyn)


def test_summarize_opt_run4(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    # Make sure info tags are handled appropriately
    atoms = bulk("Cu") * (2, 2, 1)
    atoms.info["test_dict"] = {"hi": "there", "foo": "bar"}
    atoms[0].position += [0.1, 0.1, 0.1]
    atoms.calc = EMT()
    dyn = BFGS(atoms, trajectory="test.traj")
    dyn.run()

    results = Summarize().opt(dyn)
    assert results["atoms"].info.get("test_dict", None) == {"hi": "there", "foo": "bar"}

    # test document can be jsanitized and decoded
    d = jsanitize(results, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)
    dyn.trajectory.filename = "not_a_file.traj"

    with pytest.raises(FileNotFoundError):
        Summarize().opt(dyn)


def test_vib_run1(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("N2")
    atoms.calc = EMT()
    input_atoms = deepcopy(atoms)
    vib = Vibrations(atoms)
    vib.run()

    results = VibSummarize(vib).vib()
    assert results["atoms"] == input_atoms
    assert results["molecule_metadata"]["natoms"] == len(atoms)
    assert results["parameters_vib"]["delta"] == vib.delta
    assert results["parameters_vib"]["direction"] == "central"
    assert results["parameters_vib"]["method"] == "standard"
    assert results["parameters_vib"]["ndof"] == 6
    assert results["parameters_vib"]["nfree"] == 2
    assert "nid" in results
    assert "dir_name" in results
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


def test_errors(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    with pytest.raises(
        ValueError, match="ASE Atoms object has no attached calculator."
    ):
        Summarize().run(atoms, atoms)

    initial_atoms = read(os.path.join(RUN1, "POSCAR.gz"))
    atoms = read(os.path.join(RUN1, "OUTCAR.gz"))
    atoms.calc.results = {}
    with pytest.raises(
        ValueError, match="ASE Atoms object's calculator has no results."
    ):
        Summarize().run(atoms, initial_atoms)


def test_summarize_neb(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    # Read initial and final states:
    images = read(FILE_DIR / "test_files" / "geodesic_path.xyz", index=":")
    neb = NEB(images)
    for image in images:
        image.calc = EMT()
    optimizer = BFGS(neb, trajectory="opt.traj")
    optimizer.run(fmax=0.05)

    neb_summary = Summarize().neb(
        optimizer, len(images), n_iter_return=10, trajectory=read("opt.traj", index=":")
    )
    assert neb_summary["trajectory_results"][-2]["energy"] == pytest.approx(
        1.0875330074934446, abs=1e-4
    )

    ts_atoms = neb_summary["ts_atoms"]
    ts_atoms.calc = EMT()
    assert ts_atoms.get_potential_energy() == pytest.approx(1.1603536513693768, 1e-4)
