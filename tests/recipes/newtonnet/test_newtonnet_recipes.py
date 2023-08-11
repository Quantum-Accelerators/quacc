import gzip
import os
from copy import deepcopy
from pathlib import Path
from shutil import rmtree

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule
from ase.optimize import FIRE

from quacc import SETTINGS

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None
from quacc.recipes.newtonnet.core import (
    freq_job,
    irc_job,
    quasi_irc_job,
    relax_job,
    static_job,
    ts_job,
)

try:
    import sella
except ImportError:
    sella = None

CURRENT_FILE_PATH = Path(__file__).parent.resolve()
SETTINGS.NEWTONNET_CONFIG_PATH = os.path.join(CURRENT_FILE_PATH, "config0.yml")
SETTINGS.NEWTONNET_MODEL_PATH = os.path.join(CURRENT_FILE_PATH, "best_model_state.tar")


# def startup_module():
#     """gunzip best_model_state.tar.gz with the gzip module"""
# with gzip.open(
#     os.path.join(CURRENT_FILE_PATH, "best_model_state.tar.gz"), "rb"
# ) as f_in:
#     with open(
#         os.path.join(CURRENT_FILE_PATH, "best_model_state.tar"), "wb"
#     ) as f_out:
#         f_out.write(f_in.read())


# def teardown_module():
#     """gzip back up the best_model_state.tar file with the gzip module"""
# with open(os.path.join(CURRENT_FILE_PATH, "best_model_state.tar"), "rb") as f_in:
#     with gzip.open(
#         os.path.join(CURRENT_FILE_PATH, "best_model_state.tar.gz"), "wb"
#     ) as f_out:
#         f_out.write(f_in.read())
# os.remove(os.path.join(CURRENT_FILE_PATH, "best_model_state.tar"))


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_static_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")
    output = static_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    print("energy", output["results"]["energy"])
    assert output["results"]["energy"] == pytest.approx(-9.49781055107525)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_relax_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")
    # output = relax_job(atoms, optimizer=FIRE)
    output = relax_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["results"]["energy"] == pytest.approx(-9.498925279360149)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01


xyz_data = """
18
C       0.60975077       2.99728462      -0.61541781
C      -0.59005193       2.41413160      -0.38873337
C      -0.82159168       1.01406248       0.11456567
C       0.15184655      -0.07505590      -0.39136981
C       1.54289073      -0.03399020       0.27433195
C       2.41165401       1.03626553      -0.26527396
H       0.79522243       3.35483305      -1.62935088
H       0.89415805       3.81419442       0.06095653
H      -1.49014907       3.01998005      -0.46197420
H      -0.79893963       1.00390416       1.20698553
H      -1.83871435       0.71647419      -0.15101803
H      -0.29996525      -1.03725201      -0.18869879
H       0.25992114       0.00062592      -1.46635552
H       1.41362720       0.09126494       1.35794492
H       2.01885680      -1.01246242       0.13930723
H       1.58668407       2.08346940      -0.42918377
H       3.18086711       1.34460608       0.37834213
H       2.76575603       0.84759417      -1.23589468
"""
# Define the XYZ data
xyz_data = """
8
C       0.87589219      -0.13926227       0.10160247
O      -0.07821864      -1.06257641       0.39840308
C       0.38522759      -2.67766458      -0.24094265
O       1.25339054      -2.63650946      -0.93076096
H       0.46675900       0.87142023       0.48210765
H       1.06360062      -0.03733460      -0.95079003
H       1.80266010      -0.30212484       0.61902953
H       0.92072327      -2.76503751       1.36753343
"""

# Split the XYZ data into lines
xyz_lines = xyz_data.strip().split("\n")

# Extract the number of atoms from the first line
num_atoms = int(xyz_lines[0])

# Initialize empty lists for atom symbols and positions
symbols = []
positions = []

# Parse the atom symbols and positions from the XYZ lines
for line in xyz_lines[1:]:
    parts = line.split()
    symbol = parts[0]
    x, y, z = map(float, parts[1:])
    symbols.append(symbol)
    positions.append([x, y, z])

# Create the ASE Atoms object
atoms = Atoms(symbols=symbols, positions=positions)


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_ts_job_with_default_args(tmpdir):
    tmpdir.chdir()

    # Define test inputs
    atoms = molecule("H2O")

    # Call the function
    output = ts_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "ts" in output
    assert "thermo" in output
    assert output["ts"]["results"]["energy"] == pytest.approx(-6.401303868125212)
    # assert output["ts"]["results"]["energy"] == pytest.approx(-30.935106077591225)
    # assert output["thermo"]["vib"]["results"]["imag_vib_freqs"][0] == pytest.approx(
    #     -1686.434228258355
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
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
    # traj = output["ts"]["trajectory"]

    # # Save the trajectory in XYZ format
    # write("trajectory.xyz", traj)

    # for ii, value in enumerate(traj):
    #     mlcalculator = NewtonNet(
    #         model_path=SETTINGS.NEWTONNET_MODEL_PATH.split(":"),
    #         settings_path=SETTINGS.NEWTONNET_CONFIG_PATH.split(":"),
    #     )
    #     mlcalculator.calculate(value["atoms"])
    #     energy = mlcalculator.results["energy"]
    #     energy_std = mlcalculator.results["energy_std"]
    #     outlier = mlcalculator.results["outlier"]
    #     print(
    #         f"index: {ii},"
    #         f" energy: {energy},"
    #         f" energy_std: {energy_std},"
    #         f" outlier: {outlier}"
    #     )

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "ts" in output
    assert output["ts"]["results"]["energy"] == pytest.approx(-8.897101358408396)
    assert output["thermo"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
        0.2758627558950788
    )

    # assert output["ts"]["results"]["energy"] == pytest.approx(-30.935106077591225)
    # assert output["thermo"]["vib"]["results"]["imag_vib_freqs"][0] == pytest.approx(
    #     -1686.3786461271718
    # )
    assert "thermo" in output


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_ts_job_with_custom_optimizer(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    opt_swaps = {
        "optimizer": FIRE,
    }

    # Call the function
    output = ts_job(atoms, check_convergence=False, opt_swaps=opt_swaps)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "ts" in output
    assert "thermo" in output
    assert output["ts"]["results"]["energy"] == pytest.approx(-9.498920641930049)
    # assert output["ts"]["results"]["energy"] == pytest.approx(-34.400800631142495)
    # assert output["thermo"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
    #     0.016038718562105512
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_ts_job_with_custom_optimizer_and_custom_hessian(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    opt_swaps = {
        "optimizer": FIRE,
    }

    with pytest.raises(ValueError, match="Custom hessian can only be used with Sella."):
        # Call the function
        output = ts_job(
            atoms, check_convergence=False, use_custom_hessian=True, opt_swaps=opt_swaps
        )

        # Perform assertions on the result
        assert isinstance(output, dict)
        assert "ts" in output
        assert "thermo" in output
        assert output["ts"]["results"]["energy"] == pytest.approx(-9.498920641930049)
        # assert output["ts"]["results"]["energy"] == pytest.approx(-34.400800631142495)
        # assert output["thermo"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
        #     0.016038718562105512
        # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_irc_job_with_default_args(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")

    # Call the function
    output = irc_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "thermo" in output
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.498925675881335)
    # assert output["irc"]["results"]["energy"] == pytest.approx(-34.401686)
    # assert output["thermo"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
    #     0.016611002049098378
    # )
    # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
    #     -34.401686032583896
    # )
    # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
    #     -32.658804494518634
    # )
    # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.002951301138832678
    # )
    # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -33.5387349290616
    # )
    # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
    #     1.5924339540058858
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_irc_job_with_custom_fmax(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    fmax = 0.001

    # Call the function
    output = irc_job(atoms, fmax=fmax)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "thermo" in output
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.498925675881335)
    # assert output["irc"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["thermo"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
    #     0.23465203077606747
    # )
    # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
    #     -9.41835380477994
    # )
    # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
    #     -8.713051986549818
    # )
    # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0019593027330702985
    # )
    # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -9.297218096414728
    # )
    # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
    #     0.6025061761480841
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_irc_job_with_custom_max_steps(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    max_steps = 500

    # Call the function
    output = irc_job(atoms, max_steps=max_steps)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "thermo" in output
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.498925675881335)
    # assert output["irc"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["thermo"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
    #     0.23457884629206466
    # )
    # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
    #     -9.418353073070744
    # )
    # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
    #     -8.713112216338487
    # )
    # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0019592905615021207
    # )
    # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -9.297274697250344
    # )
    # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
    #     0.6024451502493378
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_irc_job_with_custom_temperature_and_pressure(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    temperature = 500.0
    pressure = 10.0

    # Call the function
    output = irc_job(atoms, temperature=temperature, pressure=pressure)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "thermo" in output
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.498925675881335)
    # assert output["irc"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["thermo"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
    #     0.23457884629206466
    # )
    # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
    #     -9.418353073070744
    # )
    # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
    #     -8.642530002034201
    # )
    # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0019414230701050644
    # )
    # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -9.613241537086733
    # )
    # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
    #     0.6024451502493378
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
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
    assert "irc" in output
    assert "thermo" in output
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.498925675881335)
    # assert output["irc"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["thermo"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
    #     0.23457884629206466
    # )
    # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
    #     -9.418353073070744
    # )
    # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
    #     -8.713112216338487
    # )
    # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0019592905615021207
    # )
    # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -9.297274697250344
    # )
    # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
    #     0.6024451502493378
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_irc_job_with_custom_opt_swaps(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    opt_swaps = {
        "run_kwargs": {
            "direction": "reverse",
            # "direction": "forward",
            # "keep_going": True
        },
    }

    # Call the function
    output = irc_job(atoms, opt_swaps=opt_swaps)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "thermo" in output
    assert output["irc"]["results"]["energy"] == pytest.approx(-9.498925675881335)
    # assert output["irc"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["thermo"]["vib"]["results"]["vib_energies"][0] == pytest.approx(
    #     0.23457884629206466
    # )
    # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
    #     -9.418353073070744
    # )
    # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
    #     -8.713112216338487
    # )
    # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0019592905615021207
    # )
    # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -9.297274697250344
    # )
    # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
    #     0.6024451502493378
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_quasi_irc_job_with_default_args(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")

    # Call the function
    output = quasi_irc_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "opt" in output
    assert "thermo" in output
    if atoms == molecule("H2O"):
        assert output["irc"]["irc"]["results"]["energy"] == pytest.approx(
            -9.498925675881335
        )
        # assert output["irc"]["irc"]["results"]["energy"] == pytest.approx(
        #     -9.41835380477994
        # )
        # assert output["irc"]["thermo"]["vib"]["results"]["vib_energies"][
        #     0
        # ] == pytest.approx(0.23457884629206466)
        # assert output["opt"]["results"]["energy"] == pytest.approx(-9.41835380477994)
        # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
        #     -9.418353073070744
        # )
        # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
        #     -8.713112216338487
        # )
        # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
        #     0.0019592905615021207
        # )
        # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
        #     -9.297274697250344
        # )
        # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
        #     0.6024451502493378
        # )

    # assert output["irc"]["irc"]["results"]["energy"] == pytest.approx(
    #     -30.918657693685173
    # )
    # assert output["irc"]["thermo"]["vib"]["results"]["vib_energies"][
    #     0
    # ] == pytest.approx(-0.01692679363219578)
    # assert output["opt"]["results"]["energy"] == pytest.approx(-34.401864391746884)
    # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
    #     -34.401864391746884
    # )
    # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
    #     -32.65904420666086
    # )
    # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0029512963483415565
    # )
    # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -33.5389732129189
    # )
    # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
    #     1.5923690910672277
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_quasi_irc_job_with_custom_direction(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    direction = "reverse"
    # direction = "forward"

    # Call the function
    output = quasi_irc_job(atoms, direction=direction)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "opt" in output
    assert "thermo" in output
    assert output["irc"]["irc"]["results"]["energy"] == pytest.approx(
        -9.498925675881335
    )
    # assert output["irc"]["irc"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["irc"]["thermo"]["vib"]["results"]["vib_energies"][
    #     0
    # ] == pytest.approx(0.23457884629206466)
    # assert output["opt"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
    #     -9.418353073070744
    # )
    # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
    #     -8.713112216338487
    # )
    # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0019592905615021207
    # )
    # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -9.297274697250344
    # )
    # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
    #     0.6024451502493378
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_quasi_irc_job_with_custom_temperature_and_pressure(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    temperature = 500.0
    pressure = 10.0

    # Call the function
    output = quasi_irc_job(atoms, temperature=temperature, pressure=pressure)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "opt" in output
    assert "thermo" in output
    assert output["irc"]["irc"]["results"]["energy"] == pytest.approx(
        -9.498925675881335
    )
    # assert output["irc"]["irc"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["irc"]["thermo"]["vib"]["results"]["vib_energies"][
    #     0
    # ] == pytest.approx(0.23457884629206466)
    # assert output["opt"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
    #     -9.418353073070744
    # )
    # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
    #     -8.642530002034201
    # )
    # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0019414230701050644
    # )
    # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -9.613241537086733
    # )
    # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
    #     0.6024451502493378
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_quasi_irc_job_with_custom_irc_swaps(tmpdir):
    tmpdir.chdir()
    # Define test inputs
    atoms = molecule("H2O")
    irc_swaps = {
        "run_kwargs": {
            "direction": "reverse",
        },
    }

    # Call the function
    output = quasi_irc_job(atoms, irc_swaps=irc_swaps)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "opt" in output
    assert "thermo" in output
    assert output["irc"]["irc"]["results"]["energy"] == pytest.approx(
        -9.498925675881335
    )
    # assert output["irc"]["irc"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["irc"]["thermo"]["vib"]["results"]["vib_energies"][
    #     0
    # ] == pytest.approx(0.23457884629206466)
    # assert output["opt"]["results"]["energy"] == pytest.approx(-9.41835380477994)
    # assert output["thermo"]["thermo"]["results"]["energy"] == pytest.approx(
    #     -9.418353073070744
    # )
    # assert output["thermo"]["thermo"]["results"]["enthalpy"] == pytest.approx(
    #     -8.713112216338487
    # )
    # assert output["thermo"]["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0019592905615021207
    # )
    # assert output["thermo"]["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -9.297274697250344
    # )
    # assert output["thermo"]["thermo"]["results"]["zpe"] == pytest.approx(
    #     0.6024451502493378
    # )


@pytest.mark.skipif(
    NewtonNet is None or sella is None,
    reason="NewtonNet and Sella must be installed.",
)
def test_freq_job(tmpdir):
    tmpdir.chdir()
    atoms = molecule("H2O")
    output = freq_job(atoms)
    assert output["vib"]["atoms"] == molecule("H2O")
    assert len(output["vib"]["results"]["vib_freqs_raw"]) == 9
    assert len(output["vib"]["results"]["vib_freqs"]) == 3
    assert output["vib"]["results"]["vib_freqs_raw"][0] == pytest.approx(
        -5.3863100384091165e-05
    )
    # assert output["vib"]["results"]["vib_freqs_raw"][0] == pytest.approx(
    #     -292.9037211526164
    # )
    # assert output["vib"]["results"]["vib_freqs_raw"][-1] == pytest.approx(
    #     4071.5043841913816
    # )
    # assert output["vib"]["results"]["vib_freqs"][0] == pytest.approx(1863.4478302128284)
    # assert output["vib"]["results"]["vib_freqs"][-1] == pytest.approx(
    #     4071.5043841913816
    # )
    # assert output["vib"]["results"]["n_imag"] == 0
    # assert output["vib"]["results"]["imag_vib_freqs"] == []

    # assert output["thermo"]["atoms"] == atoms
    # assert output["thermo"]["symmetry"]["point_group"] == "C2v"
    # assert output["thermo"]["symmetry"]["rotation_number"] == 2
    # assert output["thermo"]["symmetry"]["linear"] is False
    # assert output["thermo"]["results"]["energy"] == -9.416288989620472
    # assert output["thermo"]["results"]["enthalpy"] == pytest.approx(-8.710428050278118)
    # assert output["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0019582532509612592
    # )
    # assert output["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -9.294281257052218
    # )

    atoms = molecule("CH3")
    initial_atoms = deepcopy(atoms)
    output = freq_job(atoms, temperature=1000, pressure=20)
    assert output["vib"]["atoms"] == initial_atoms
    assert len(output["vib"]["results"]["vib_freqs_raw"]) == 12
    assert len(output["vib"]["results"]["vib_freqs"]) == 6
    assert output["vib"]["results"]["vib_energies_raw"][0] == pytest.approx(
        -0.08792854890509715
    )
    # assert output["vib"]["results"]["vib_energies_raw"][0] == pytest.approx(
    #     -1.5636131383185185e-09
    # )
    # assert output["vib"]["results"]["vib_energies_raw"][-1] == pytest.approx(
    #     0.40918939605228777
    # )
    # assert output["vib"]["results"]["vib_energies"][0] == pytest.approx(
    #     0.09310183598706114
    # )
    # assert output["vib"]["results"]["vib_energies"][-1] == pytest.approx(
    #     0.40918939605228777
    # )
    # assert output["vib"]["results"]["n_imag"] == 0
    # assert output["vib"]["results"]["imag_vib_freqs"] == []
    # assert output["thermo"]["atoms"] == initial_atoms
    # assert output["thermo"]["atoms"] == initial_atoms
