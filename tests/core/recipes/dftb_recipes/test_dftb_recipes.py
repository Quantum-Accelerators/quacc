from __future__ import annotations

from shutil import which

import pytest

DFTBPLUS_EXISTS = bool(which("dftb+"))

pytestmark = pytest.mark.skipif(not DFTBPLUS_EXISTS, reason="Needs DFTB+")
import os
from logging import getLogger

import numpy as np
from ase.build import bulk, molecule

from quacc import JobFailure, change_settings
from quacc.recipes.dftb.core import relax_job, static_job

LOGGER = getLogger(__name__)
LOGGER.propagate = True


def test_static_job_water(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    atoms.info = {"test": "hello"}
    output = static_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert output["parameters"]["Hamiltonian_MaxSccIterations"] == 200
    assert output["results"]["energy"] == pytest.approx(-137.9677759924738)
    assert (
        np.array_equal(output["atoms"].get_positions(), atoms.get_positions()) is True
    )
    assert output["atoms"].info["test"] == "hello"
    assert output["atoms"].info.get("_id")


def test_static_job_cu_supercell(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu") * (3, 3, 3)
    output = static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert output["parameters"]["Hamiltonian_MaxSccIterations"] == 200
    assert (
        output["parameters"]["Hamiltonian_KPointsAndWeights_"].strip()
        == "SupercellFolding"
    )
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty000"] == "1 0 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty001"] == "0 1 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty002"] == "0 0 1"
    assert output["results"]["energy"] == pytest.approx(-2885.417379886678)
    assert (
        np.array_equal(output["atoms"].get_positions(), atoms.get_positions()) is True
    )
    assert np.array_equal(output["atoms"].cell.array, atoms.cell.array) is True


def test_static_job_cu_kpts(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    output = static_job(atoms, kpts=(3, 3, 3))
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert output["parameters"]["Hamiltonian_MaxSccIterations"] == 200
    assert (
        output["parameters"]["Hamiltonian_KPointsAndWeights_"].strip()
        == "SupercellFolding"
    )
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty000"] == "3 0 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty001"] == "0 3 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty002"] == "0 0 3"
    assert output["results"]["energy"] == pytest.approx(-106.86647125470942)
    assert (
        np.array_equal(output["atoms"].get_positions(), atoms.get_positions()) is True
    )
    assert np.array_equal(output["atoms"].cell.array, atoms.cell.array) is True


def test_static_errors(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = molecule("H2O")
    with pytest.raises(JobFailure, match="Calculation failed!") as err:
        static_job(atoms, Hamiltonian_MaxSccIterations=1)
    with pytest.raises(RuntimeError, match="failed with command"):
        raise err.value.parent_error


def test_relax_job_water(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")

    output = relax_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert output["parameters"]["Hamiltonian_MaxSccIterations"] == 200
    assert output["results"]["energy"] == pytest.approx(-137.97654214864497)
    assert (
        np.array_equal(output["atoms"].get_positions(), atoms.get_positions()) is False
    )
    assert np.array_equal(output["atoms"].cell.array, atoms.cell.array) is True


def test_relax_job_cu_supercell(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1

    output = relax_job(atoms, kpts=(3, 3, 3), Hamiltonian_MaxSccIterations=100)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert output["parameters"]["Hamiltonian_MaxSccIterations"] == 100
    assert (
        output["parameters"]["Hamiltonian_KPointsAndWeights_"].strip()
        == "SupercellFolding"
    )
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty000"] == "3 0 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty001"] == "0 3 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty002"] == "0 0 3"
    assert output["parameters"]["Driver_"] == "GeometryOptimization"
    assert output["parameters"]["Driver_LatticeOpt"] == "No"
    assert (
        np.array_equal(output["atoms"].get_positions(), atoms.get_positions()) is False
    )
    assert np.array_equal(output["atoms"].cell.array, atoms.cell.array) is True


def test_relax_job_cu_supercell_cell_relax(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    output = relax_job(atoms, method="GFN1-xTB", kpts=(3, 3, 3), relax_cell=True)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN1-xTB"
    assert output["parameters"]["Hamiltonian_MaxSccIterations"] == 200
    assert (
        output["parameters"]["Hamiltonian_KPointsAndWeights_"].strip()
        == "SupercellFolding"
    )
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty000"] == "3 0 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty001"] == "0 3 0"
    assert output["parameters"]["Hamiltonian_KPointsAndWeights_empty002"] == "0 0 3"
    assert output["parameters"]["Driver_"] == "GeometryOptimization"
    assert output["parameters"]["Driver_LatticeOpt"] == "Yes"
    assert (
        np.array_equal(output["atoms"].get_positions(), atoms.get_positions()) is False
    )
    assert np.array_equal(output["atoms"].cell.array, atoms.cell.array) is False


def test_relax_job_cu_supercell_errors(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.5
    with pytest.raises(JobFailure, match="Calculation failed!") as err:
        relax_job(atoms, kpts=(3, 3, 3), MaxSteps=1, Hamiltonian_MaxSccIterations=100)
    with pytest.raises(RuntimeError, match="failed with command"):
        raise err.value.parent_error


@pytest.mark.skipif(os.name == "nt", reason="symlinking not possible on Windows")
def test_child_errors(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    with (
        caplog.at_level(LOGGER.INFO),
        change_settings({"SCRATCH_DIR": tmp_path / "scratch"}),
    ):
        with pytest.raises(JobFailure, match="Calculation failed!") as err:
            static_job(atoms)
        with pytest.raises(RuntimeError, match="failed with command"):
            raise err.value.parent_error
        assert "Calculation failed" in caplog.text
        assert "failed-quacc-" in " ".join(os.listdir(tmp_path / "scratch"))


@pytest.mark.skipif(os.name == "nt", reason="symlinking not possible on Windows")
def test_child_errors2(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    atoms = bulk("Cu")
    with (
        caplog.at_level(LOGGER.INFO),
        change_settings({"SCRATCH_DIR": tmp_path / "scratch"}),
    ):
        with pytest.raises(JobFailure, match="Calculation failed!") as err:
            relax_job(atoms)
        with pytest.raises(RuntimeError, match="failed with command"):
            raise err.value.parent_error
        assert "Calculation failed" in caplog.text
        assert "failed-quacc-" in " ".join(os.listdir(tmp_path / "scratch"))
