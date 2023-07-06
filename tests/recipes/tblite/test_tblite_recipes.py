import os
from copy import deepcopy
from shutil import rmtree

import numpy as np
import pytest
from ase.build import molecule

from quacc.recipes.tblite.core import freq_job, relax_job, static_job

try:
    from tblite.ase import TBLite
except ImportError:
    TBLite = None

try:
    import sella
except ImportError:
    sella = None


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
    TBLite is None,
    reason="tblite must be installed.",
)
def teststatic_job():
    atoms = molecule("H2O")
    output = static_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.96777594361672)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    output = static_job(atoms, method="GFN1-xTB")
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["results"]["energy"] == pytest.approx(-156.96750578831137)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


@pytest.mark.skipif(
    TBLite is None,
    reason="tblite must be installed.",
)
def test_relax_Job():
    atoms = molecule("H2O")
    output = relax_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.97654191396492)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01


@pytest.mark.skipif(
    TBLite is None,
    reason="tblite must be installed.",
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
    assert len(output["thermo"]["parameters_thermo"]["vib_freqs"]) == 3
    assert output["vib"]["results"]["vib_freqs"][0] == pytest.approx(1586.623114694335)
    assert output["thermo"]["parameters_thermo"]["vib_freqs"][-1] == pytest.approx(
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
    assert len(output["thermo"]["parameters_thermo"]["vib_freqs"]) == 0
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
    assert output["thermo"]["parameters_thermo"]["temperature"] == 1000.0
    assert output["thermo"]["parameters_thermo"]["pressure"] == 20.0
    assert output["thermo"]["parameters_thermo"]["sigma"] == 6
    assert output["thermo"]["parameters_thermo"]["spin_multiplicity"] == 2
    assert output["thermo"]["symmetry"]["linear"] is False
    assert output["thermo"]["symmetry"]["rotation_number"] == 6
    assert len(output["thermo"]["parameters_thermo"]["vib_freqs"]) == 6
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
