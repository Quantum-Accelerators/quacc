from shutil import which

import numpy as np
import pytest
from ase.build import bulk, molecule

from quacc import SETTINGS
from quacc.recipes.dftb.core import relax_job, static_job

DFTBPLUS_EXISTS = bool(which("dftb+"))


@pytest.mark.skipif(
    DFTBPLUS_EXISTS is False,
    reason="DFTB+ must be installed. Try conda install -c conda-forge dftbplus",
)
def test_static_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")

    output = static_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.9677759924738)
    assert (
        np.array_equal(output["atoms"].get_positions(), atoms.get_positions()) is True
    )

    atoms = bulk("Cu")

    output = static_job(atoms, kpts=(3, 3, 3))
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
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

    with pytest.raises(ValueError):
        atoms = molecule("H2O")
        output = static_job(atoms, calc_swaps={"Hamiltonian_MaxSccIterations": 1})


@pytest.mark.skipif(
    DFTBPLUS_EXISTS is False,
    reason="DFTB+ must be installed. Try conda install -c conda-forge dftbplus",
)
def test_relax_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")

    output = relax_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.97654214864497)
    assert (
        np.array_equal(output["atoms"].get_positions(), atoms.get_positions()) is False
    )
    assert np.array_equal(output["atoms"].cell.array, atoms.cell.array) is True

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1

    output = relax_job(atoms, kpts=(3, 3, 3))
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN2-xTB"
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

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    output = relax_job(atoms, method="GFN1-xTB", kpts=(3, 3, 3), lattice_opt=True)
    assert output["nsites"] == len(atoms)
    assert output["parameters"]["Hamiltonian_"] == "xTB"
    assert output["parameters"]["Hamiltonian_Method"] == "GFN1-xTB"
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

    with pytest.raises(ValueError):
        atoms = bulk("Cu") * (2, 1, 1)
        atoms[0].position += 0.5
        relax_job(atoms, kpts=(3, 3, 3), calc_swaps={"MaxSteps": 1})


@pytest.mark.skipif(
    DFTBPLUS_EXISTS is False,
    reason="DFTB+ must be installed. Try conda install -c conda-forge dftbplus",
)
def test_unique_workdir(tmpdir):
    SETTINGS.CREATE_UNIQUE_WORKDIR = True
    test_static_job(tmpdir)
    test_relax_job(tmpdir)
    SETTINGS.CREATE_UNIQUE_WORKDIR = DEFAULT_SETTINGS.CREATE_UNIQUE_WORKDIR
