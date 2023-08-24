import numpy as np
import pytest
from ase.build import molecule

from quacc import SETTINGS
from quacc.recipes.lj.core import freq_job, relax_job, static_job


@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_static_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")

    output = static_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["epsilon"] == 1.0
    assert output["parameters"]["sigma"] == 1.0
    assert output["parameters"]["rc"] == 3
    assert output["parameters"]["ro"] == 0.66 * 3
    assert output["results"]["energy"] == pytest.approx(1.772068860679255)

    output = static_job(atoms, calc_swaps={"epsilon": 2.0, "rc": 0.5})
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["epsilon"] == 2.0
    assert output["parameters"]["sigma"] == 1.0
    assert output["parameters"]["rc"] == 0.5
    assert output["parameters"]["ro"] == 0.66 * 0.5
    assert output["results"]["energy"] == pytest.approx(0.0)


@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_relax_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")

    output = relax_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["epsilon"] == 1.0
    assert output["parameters"]["sigma"] == 1.0
    assert output["parameters"]["rc"] == 3
    assert output["parameters"]["ro"] == 0.66 * 3
    assert output["results"]["energy"] == pytest.approx(-2.983561029599189)
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01

    output = relax_job(
        atoms, opt_swaps={"fmax": 0.03}, calc_swaps={"epsilon": 2.0, "rc": 0.5}
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["epsilon"] == 2.0
    assert output["parameters"]["sigma"] == 1.0
    assert output["parameters"]["rc"] == 0.5
    assert output["parameters"]["ro"] == 0.66 * 0.5
    assert output["results"]["energy"] == pytest.approx(0.0)
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.03


@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_freq_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")

    output = freq_job(relax_job(atoms))
    assert output["vib"]["natoms"] == len(atoms)
    assert output["thermo"]["natoms"] == len(atoms)
    assert output["vib"]["parameters"]["epsilon"] == 1.0
    assert output["vib"]["parameters"]["sigma"] == 1.0
    assert output["vib"]["parameters"]["rc"] == 3
    assert output["vib"]["parameters"]["ro"] == 0.66 * 3
    assert len(output["vib"]["results"]["vib_freqs_raw"]) == 3 * len(atoms)
    assert len(output["vib"]["results"]["vib_freqs"]) == 3 * len(atoms) - 6
    assert len(output["thermo"]["parameters_thermo"]["vib_freqs"]) == 3 * len(atoms) - 6
    assert output["thermo"]["parameters_thermo"]["n_imag"] == 0
