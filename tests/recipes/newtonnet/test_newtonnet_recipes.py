import os
from pathlib import Path

import numpy as np
import pytest
from ase.build import molecule
from ase.optimize import FIRE

from quacc import SETTINGS
from quacc.recipes.newtonnet.core import freq_job, relax_job, static_job
from quacc.recipes.newtonnet.ts import irc_job, quasi_irc_job, ts_job

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None

try:
    import sella
except ImportError:
    sella = None

CURRENT_FILE_PATH = Path(__file__).parent.resolve()
SETTINGS.NEWTONNET_CONFIG_PATH = os.path.join(CURRENT_FILE_PATH, "config0.yml")
SETTINGS.NEWTONNET_MODEL_PATH = os.path.join(CURRENT_FILE_PATH, "best_model_state.tar")


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_static_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")
    output = static_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["results"]["energy"] == pytest.approx(-9.515200426406743)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_relax_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")
    # output = relax_job(atoms, optimizer=FIRE)
    output = relax_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["results"]["energy"] == pytest.approx(-9.517354818364769)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_freq_job(tmpdir):
    tmpdir.chdir()
    atoms = molecule("H2O")
    output = freq_job(atoms)
    assert output["atoms"] == molecule("H2O")
    assert len(output["results"]["vib_freqs_raw"]) == 9
    assert len(output["results"]["vib_freqs"]) == 3
    assert output["results"]["vib_freqs_raw"][-1] == pytest.approx(4090.37777396351)
    assert output["results"]["vib_freqs"][0] == pytest.approx(1814.0941260498644)
    assert output["results"]["vib_freqs"][-1] == pytest.approx(4090.37777396351)
    assert output["results"]["n_imag"] == 0
    assert output["results"]["imag_vib_freqs"] == []

    assert output["thermo"]["atoms"] == atoms
    assert output["thermo"]["symmetry"]["point_group"] == "C2v"
    assert output["thermo"]["symmetry"]["rotation_number"] == 2
    assert output["thermo"]["symmetry"]["linear"] is False
    assert output["thermo"]["results"]["energy"] == pytest.approx(-9.515200426406743)
    assert output["thermo"]["results"]["enthalpy"] == pytest.approx(-8.807932688921495)
    assert output["thermo"]["results"]["entropy"] == pytest.approx(
        0.0019582788098945945
    )
    assert output["thermo"]["results"]["gibbs_energy"] == pytest.approx(
        -9.391793516091568
    )

    atoms = molecule("CH3")
    output = freq_job(atoms, temperature=1000, pressure=20)
    assert output["atoms"] == molecule("CH3")
    assert len(output["results"]["vib_freqs_raw"]) == 12
    assert len(output["results"]["vib_freqs"]) == 6
    assert output["results"]["vib_energies_raw"][0] == pytest.approx(
        -0.09441402482739979
    )
    assert output["results"]["vib_energies_raw"][-1] == pytest.approx(
        0.3925829460532815
    )
    assert output["results"]["vib_energies"][0] == pytest.approx(-0.09441402482739979)
    assert output["results"]["vib_energies"][-1] == pytest.approx(0.3925829460532815)
    assert output["results"]["n_imag"] == 1
    assert output["results"]["imag_vib_freqs"] == pytest.approx([-761.5004719152678])
    assert output["thermo"]["atoms"] == molecule("CH3")


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_ts_job_with_default_args(tmpdir):
    tmpdir.chdir()

    # Define test inputs
    atoms = molecule("H2O")

    # Call the function
    output = ts_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)

    assert "freq" in output
    assert "thermo" in output["freq"]
    assert output["results"]["energy"] == pytest.approx(-6.796914263061945)
    assert output["freq"]["vib"]["results"]["imag_vib_freqs"][0] == pytest.approx(
        -2426.7398321816004
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_ts_job_with_custom_hessian(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    use_custom_hessian = True
    opt_swaps = {"max_steps": 4}
    # Call the function
    output = ts_job(
        atoms,
        use_custom_hessian=use_custom_hessian,
        check_convergence=False,
        opt_swaps=opt_swaps,
    )

    # Perform assertions on the result
    assert isinstance(output, dict)

    assert output["results"]["energy"] == pytest.approx(-8.855604432470276)
    assert output["freq"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
        0.2256022513686731
    )
    assert "thermo" in output["freq"]


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_ts_job_with_custom_optimizer(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    opt_swaps = {"optimizer": FIRE}

    # Call the function
    output = ts_job(atoms, check_convergence=False, opt_swaps=opt_swaps)

    # Perform assertions on the result
    assert isinstance(output, dict)

    assert "thermo" in output["freq"]
    assert output["results"]["energy"] == pytest.approx(-9.51735515322368)
    assert output["freq"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
        0.22679888726664774
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_ts_job_with_custom_optimizer_and_custom_hessian(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    opt_swaps = {"optimizer": FIRE}

    with pytest.raises(ValueError, match="Custom hessian can only be used with Sella."):
        # Call the function
        ts_job(
            atoms, check_convergence=False, use_custom_hessian=True, opt_swaps=opt_swaps
        )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_irc_job_with_default_args(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")

    # Call the function
    output = irc_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "vib" in output["freq"]
    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq"]["thermo"]["results"]["energy"] == pytest.approx(
        -9.517354091813969
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_irc_job_with_custom_fmax(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    fmax = 0.001

    # Call the function
    output = irc_job(atoms, opt_swaps={"fmax": fmax})

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "vib" in output["freq"]
    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq"]["thermo"]["results"]["energy"] == pytest.approx(
        -9.517354091813969
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_irc_job_with_custom_max_steps(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    max_steps = 500

    # Call the function
    output = irc_job(atoms, opt_swaps={"max_steps": max_steps})

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "vib" in output["freq"]
    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq"]["thermo"]["results"]["energy"] == pytest.approx(
        -9.517354091813969
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_irc_job_with_custom_temperature_and_pressure(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    temperature = 500.0
    pressure = 10.0

    # Call the function
    output = irc_job(
        atoms, freq_job_kwargs={"temperature": temperature, "pressure": pressure}
    )

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "vib" in output["freq"]
    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq"]["thermo"]["results"]["energy"] == pytest.approx(
        -9.517354091813969
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_irc_job_with_check_convergence(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    check_convergence = True

    # Call the function
    output = irc_job(atoms, check_convergence=check_convergence)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "vib" in output["freq"]
    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq"]["thermo"]["results"]["energy"] == pytest.approx(
        -9.517354091813969
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_irc_job_with_custom_opt_swaps(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    opt_swaps = {"run_kwargs": {"direction": "reverse"}}

    # Call the function
    output = irc_job(atoms, opt_swaps=opt_swaps)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "vib" in output["freq"]
    assert output["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["freq"]["thermo"]["results"]["energy"] == pytest.approx(
        -9.517354965639784
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_quasi_irc_job_with_default_args(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")

    # Call the function
    output = quasi_irc_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "vib" in output["freq"]
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq"]["thermo"]["results"]["energy"] == pytest.approx(
        -9.517354091813969
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_quasi_irc_job_with_custom_direction(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    direction = "reverse"

    # Call the function
    output = quasi_irc_job(atoms, direction=direction)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "vib" in output["freq"]
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["freq"]["thermo"]["results"]["energy"] == pytest.approx(
        -9.517354965639784
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_quasi_irc_job_with_custom_temperature_and_pressure(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    temperature = 500.0
    pressure = 10.0

    # Call the function
    output = quasi_irc_job(
        atoms, freq_job_kwargs={"temperature": temperature, "pressure": pressure}
    )

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "vib" in output["freq"]
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["results"]["energy"] == pytest.approx(-9.517354091813969)
    assert output["freq"]["thermo"]["results"]["energy"] == pytest.approx(
        -9.517354091813969
    )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {None, "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_quasi_irc_job_with_custom_irc_swaps(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    irc_swaps = {
        "run_kwargs": {"direction": "reverse"},
    }

    # Call the function
    output = quasi_irc_job(atoms, irc_swaps=irc_swaps)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "vib" in output["freq"]
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["results"]["energy"] == pytest.approx(-9.517354965639784)
    assert output["freq"]["thermo"]["results"]["energy"] == pytest.approx(
        -9.517354965639784
    )
