import os
from copy import deepcopy
from shutil import rmtree

import numpy as np
import pytest
from ase.build import molecule
from ase.optimize import FIRE
from ase import Atoms


try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None
from quacc.recipes.newtonnet.core import ts_job, irc_job, quasi_irc_job, freq_job, relax_job, static_job


def teardown_module():
    for f in os.listdir("."):
        if (
            ".log" in f
            or ".pckl" in f
            or ".traj" in f
            or "gfnff_topo" in f
            or ".gz" in f
        ):
            os.remove(f)
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir" or f == "vib":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


@pytest.mark.skipif(
    NewtonNet is None,
    reason="NewtonNet must be installed.",
)
def test_static_Job():
    atoms = molecule("H2O")
    output = static_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    print('energy' ,output['results']['energy'])
    assert output["results"]["energy"] == pytest.approx(-9.416288989620472)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


@pytest.mark.skipif(
    NewtonNet is None,
    reason="NewtonNet must be installed.",
)
def test_relax_Job():
    atoms = molecule("H2O")
    #output = relax_job(atoms, optimizer=FIRE)
    output = relax_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["results"]["energy"] == pytest.approx(-9.418353066770434)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01


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
xyz_lines = xyz_data.strip().split('\n')

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


def test_ts_job_with_default_args():
    # Define test inputs
    #atoms = molecule("H2O")

    # Call the function
    output = ts_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "ts" in output
    assert "thermo" in output
    assert output['ts']['results']['energy'] == pytest.approx(-30.935106077591225)
    assert output['thermo']['vib']['results']['imag_vib_freqs'][0] == pytest.approx(-1686.434228258355)

def test_ts_job_with_custom_hessian():
    # Define test inputs
    #atoms = molecule("H2O")
    use_custom_hessian = True

    # Call the function
    output = ts_job(atoms, use_custom_hessian=use_custom_hessian)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "ts" in output
    assert output['ts']['results']['energy'] == pytest.approx(-30.935106077591225)
    assert output['thermo']['vib']['results']['imag_vib_freqs'][0] == pytest.approx(-1686.3786461271718)
    assert "thermo" in output


def test_ts_job_with_custom_optimizer():
    # Define test inputs
    # atoms = molecule("H2O")
    opt_swaps = {
        "optimizer": FIRE,
    }

    # Call the function
    output = ts_job(atoms, opt_swaps=opt_swaps)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "ts" in output
    assert "thermo" in output
    assert output['ts']['results']['energy'] == pytest.approx(-34.400800631142495)
    assert output['thermo']['vib']['results']['vib_energies'][0] == pytest.approx(0.016038718562105512)


def test_irc_job_with_default_args():
    # Define test inputs
    # atoms = molecule("H2O")

    # Call the function
    output = irc_job(atoms)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "thermo" in output
    assert output['irc']['results']['energy'] == pytest.approx(-34.401686)
    assert output['thermo']['vib']['results']['vib_energies'][0] == pytest.approx(0.016611002049098378)
    assert output['thermo']['thermo']['results']['energy'] == pytest.approx(-34.401686032583896)
    assert output['thermo']['thermo']['results']['enthalpy'] == pytest.approx(-32.658804494518634)
    assert output['thermo']['thermo']['results']['entropy'] == pytest.approx(0.002951301138832678)
    assert output['thermo']['thermo']['results']['gibbs_energy'] == pytest.approx(-33.5387349290616)
    assert output['thermo']['thermo']['results']['zpe'] == pytest.approx(1.5924339540058858)


def test_irc_job_with_custom_fmax():
    # Define test inputs
    atoms = molecule("H2O")
    fmax = 0.001

    # Call the function
    output = irc_job(atoms, fmax=fmax)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "thermo" in output
    assert output['irc']['results']['energy'] == pytest.approx(-9.41835380477994)
    assert output['thermo']['vib']['results']['vib_energies'][0] == pytest.approx(0.23465203077606747)
    assert output['thermo']['thermo']['results']['energy'] == pytest.approx(-9.41835380477994)
    assert output['thermo']['thermo']['results']['enthalpy'] == pytest.approx(-8.713051986549818)
    assert output['thermo']['thermo']['results']['entropy'] == pytest.approx(0.0019593027330702985)
    assert output['thermo']['thermo']['results']['gibbs_energy'] == pytest.approx(-9.297218096414728)
    assert output['thermo']['thermo']['results']['zpe'] == pytest.approx(0.6025061761480841)


def test_irc_job_with_custom_max_steps():
    # Define test inputs
    atoms = molecule("H2O")
    max_steps = 500

    # Call the function
    output = irc_job(atoms, max_steps=max_steps)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "thermo" in output
    assert output['irc']['results']['energy'] == pytest.approx(-9.41835380477994)
    assert output['thermo']['vib']['results']['vib_energies'][0] == pytest.approx(0.23457884629206466)
    assert output['thermo']['thermo']['results']['energy'] == pytest.approx(-9.418353073070744)
    assert output['thermo']['thermo']['results']['enthalpy'] == pytest.approx(-8.713112216338487)
    assert output['thermo']['thermo']['results']['entropy'] == pytest.approx(0.0019592905615021207)
    assert output['thermo']['thermo']['results']['gibbs_energy'] == pytest.approx(-9.297274697250344)
    assert output['thermo']['thermo']['results']['zpe'] == pytest.approx(0.6024451502493378)


def test_irc_job_with_custom_temperature_and_pressure():
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
    assert output['irc']['results']['energy'] == pytest.approx(-9.41835380477994)
    assert output['thermo']['vib']['results']['vib_energies'][0] == pytest.approx(0.23457884629206466)
    assert output['thermo']['thermo']['results']['energy'] == pytest.approx(-9.418353073070744)
    assert output['thermo']['thermo']['results']['enthalpy'] == pytest.approx(-8.642530002034201)
    assert output['thermo']['thermo']['results']['entropy'] == pytest.approx(0.0019414230701050644)
    assert output['thermo']['thermo']['results']['gibbs_energy'] == pytest.approx(-9.613241537086733)
    assert output['thermo']['thermo']['results']['zpe'] == pytest.approx(0.6024451502493378)


def test_irc_job_with_check_convergence():
    # Define test inputs
    atoms = molecule("H2O")
    check_convergence = True

    # Call the function
    output = irc_job(atoms, check_convergence=check_convergence)

    # Perform assertions on the result
    assert isinstance(output, dict)
    assert "irc" in output
    assert "thermo" in output
    assert output['irc']['results']['energy'] == pytest.approx(-9.41835380477994)
    assert output['thermo']['vib']['results']['vib_energies'][0] == pytest.approx(0.23457884629206466)
    assert output['thermo']['thermo']['results']['energy'] == pytest.approx(-9.418353073070744)
    assert output['thermo']['thermo']['results']['enthalpy'] == pytest.approx(-8.713112216338487)
    assert output['thermo']['thermo']['results']['entropy'] == pytest.approx(0.0019592905615021207)
    assert output['thermo']['thermo']['results']['gibbs_energy'] == pytest.approx(-9.297274697250344)
    assert output['thermo']['thermo']['results']['zpe'] == pytest.approx(0.6024451502493378)

'''
def test_irc_job_with_custom_newtonnet_kwargs():
    # Define test inputs
    atoms = molecule("H2O")
    newtonnet_kwargs = {
        "model_path": "/path/to/custom/model",
        "settings_path": "/path
'''


'''
@pytest.mark.skipif(
    NewtonNet is None,
    reason="NewtonNet must be installed.",
)
def test_freq_job():
    atoms = molecule("H2O")
    output = freq_job(atoms)
    assert output["vib"]["atoms"] == molecule("H2O")
    assert len(output["vib"]["results"]["vib_freqs_raw"]) == 9
    assert len(output["vib"]["results"]["vib_freqs"]) == 3
    assert output["vib"]["results"]["vib_freqs_raw"][0] == pytest.approx(
        -0.10864429415434408
    )
    assert output["vib"]["results"]["vib_freqs_raw"][-1] == pytest.approx(
        3526.9940431752034
    )
    assert output["vib"]["results"]["vib_freqs"][0] == pytest.approx(1586.623114694335)
    assert output["vib"]["results"]["vib_freqs"][-1] == pytest.approx(
        3526.9940431752034
    )
    assert output["vib"]["results"]["n_imag"] == 0
    assert output["vib"]["results"]["imag_vib_freqs"] == []

    assert output["thermo"]["atoms"] == atoms
    assert output["thermo"]["symmetry"]["point_group"] == "C2v"
    assert output["thermo"]["symmetry"]["rotation_number"] == 2
    assert output["thermo"]["symmetry"]["linear"] is False
    assert len(output["thermo"]["results"]["vib_freqs"]) == 3
    assert output["vib"]["results"]["vib_freqs"][0] == pytest.approx(1586.623114694335)
    assert output["thermo"]["results"]["vib_freqs"][-1] == pytest.approx(
        3526.9940431752034
    )
    assert output["thermo"]["results"]["energy"] == 0.0
    assert output["thermo"]["results"]["enthalpy"] == pytest.approx(0.6375973622705722)
    assert output["thermo"]["results"]["entropy"] == pytest.approx(
        0.0019584992229988523
    )
    assert output["thermo"]["results"]["gibbs_energy"] == pytest.approx(
        0.05367081893346437
    )

    atoms = molecule("H")
    atoms.set_initial_magnetic_moments([0.0])
    initial_atoms = deepcopy(atoms)
    output = freq_job(atoms, energy=-1.0)
    assert output["vib"]["atoms"] == initial_atoms
    assert len(output["vib"]["results"]["vib_freqs_raw"]) == 3
    assert len(output["vib"]["results"]["vib_freqs"]) == 0
    assert output["vib"]["results"]["vib_freqs_raw"][0] == 0
    assert output["vib"]["results"]["vib_freqs_raw"][-1] == 0
    assert output["vib"]["results"]["vib_freqs"] == []
    assert output["vib"]["results"]["n_imag"] == 0
    assert output["vib"]["results"]["imag_vib_freqs"] == []

    assert output["thermo"]["atoms"] == initial_atoms
    assert output["thermo"]["symmetry"]["linear"] is False
    assert output["thermo"]["symmetry"]["rotation_number"] == np.inf
    assert len(output["thermo"]["results"]["vib_freqs"]) == 0
    assert output["thermo"]["results"]["energy"] == -1.0
    assert output["thermo"]["results"]["enthalpy"] == pytest.approx(-0.9357685739989672)
    assert output["thermo"]["results"]["entropy"] == pytest.approx(
        0.0011292352752446438
    )
    assert output["thermo"]["results"]["gibbs_energy"] == pytest.approx(
        -1.2724500713131577
    )

    atoms = molecule("CH3")
    initial_atoms = deepcopy(atoms)
    output = freq_job(atoms, energy=-10.0, temperature=1000, pressure=20)
    assert output["vib"]["atoms"] == initial_atoms
    assert len(output["vib"]["results"]["vib_freqs_raw"]) == 12
    assert len(output["vib"]["results"]["vib_freqs"]) == 6
    assert output["vib"]["results"]["vib_energies_raw"][0] == pytest.approx(
        -9.551076713062095e-06
    )
    assert output["vib"]["results"]["vib_energies_raw"][-1] == pytest.approx(
        0.3880868821616259
    )
    assert output["vib"]["results"]["vib_energies"][0] == pytest.approx(
        0.0713506770137291
    )
    assert output["vib"]["results"]["vib_energies"][-1] == pytest.approx(
        0.3880868821616259
    )
    assert output["vib"]["results"]["n_imag"] == 0
    assert output["vib"]["results"]["imag_vib_freqs"] == []

    assert output["thermo"]["atoms"] == initial_atoms
    assert output["thermo"]["parameters"]["temperature"] == 1000.0
    assert output["thermo"]["parameters"]["pressure"] == 20.0
    assert output["thermo"]["parameters"]["sigma"] == 6
    assert output["thermo"]["parameters"]["spin_multiplicity"] == 2
    assert output["thermo"]["symmetry"]["linear"] is False
    assert output["thermo"]["symmetry"]["rotation_number"] == 6
    assert len(output["thermo"]["results"]["vib_freqs"]) == 6
    assert output["thermo"]["results"]["energy"] == -10.0
    assert output["thermo"]["results"]["enthalpy"] == pytest.approx(-8.749341973959462)
    assert output["thermo"]["results"]["entropy"] == pytest.approx(
        0.0023506788982171896
    )
    assert output["thermo"]["results"]["gibbs_energy"] == pytest.approx(
        -11.100020872176652
    )
    assert "nid" in output["thermo"]
    assert "dir_name" in output["thermo"]
    assert "nid" in output["vib"]
    assert "dir_name" in output["vib"]
'''
