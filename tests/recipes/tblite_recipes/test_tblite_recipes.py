import pytest

from quacc import SETTINGS

pytest.importorskip("tblite.ase")


DEFAULT_SETTINGS = SETTINGS.copy()


def setup_module():
    SETTINGS.WORKFLOW_ENGINE = "local"


def teardown_module():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE


def test_static_job(tmpdir):
    import numpy as np
    from ase.build import molecule

    from quacc.recipes.tblite.core import static_job

    tmpdir.chdir()

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


def test_relax_job(tmpdir):
    import numpy as np
    from ase.build import molecule

    from quacc.recipes.tblite.core import relax_job

    tmpdir.chdir()

    atoms = molecule("H2O")
    output = relax_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.97654191396492)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())
    assert np.max(np.linalg.norm(output["results"]["forces"], axis=1)) < 0.01


def test_freq_job(tmpdir):
    from copy import deepcopy

    import numpy as np
    from ase.build import molecule

    from quacc.recipes.tblite.core import freq_job

    tmpdir.chdir()

    atoms = molecule("H2O")
    output = freq_job(atoms)
    assert output["atoms"] == molecule("H2O")
    assert len(output["results"]["vib_freqs_raw"]) == 9
    assert len(output["results"]["vib_freqs"]) == 3
    assert output["results"]["vib_freqs_raw"][0] == pytest.approx(-0.10864429415434408)
    assert output["results"]["vib_freqs_raw"][-1] == pytest.approx(3526.9940431752034)
    assert output["results"]["vib_freqs"][0] == pytest.approx(1586.623114694335)
    assert output["results"]["vib_freqs"][-1] == pytest.approx(3526.9940431752034)
    assert output["results"]["n_imag"] == 0
    assert output["results"]["imag_vib_freqs"] == []

    assert output["symmetry"]["point_group"] == "C2v"
    assert output["symmetry"]["rotation_number"] == 2
    assert output["symmetry"]["linear"] is False
    assert len(output["parameters_thermo"]["vib_freqs"]) == 3
    assert output["results"]["vib_freqs"][0] == pytest.approx(1586.623114694335)
    assert output["parameters_thermo"]["vib_freqs"][-1] == pytest.approx(
        3526.9940431752034
    )
    assert output["results"]["energy"] == 0.0
    assert output["results"]["enthalpy"] == pytest.approx(0.6375973622705722)
    assert output["results"]["entropy"] == pytest.approx(0.0019584992229988523)
    assert output["results"]["gibbs_energy"] == pytest.approx(0.05367081893346437)

    atoms = molecule("H")
    atoms.set_initial_magnetic_moments([0.0])
    initial_atoms = deepcopy(atoms)
    output = freq_job(atoms, energy=-1.0)
    assert output["atoms"] == initial_atoms
    assert len(output["results"]["vib_freqs_raw"]) == 3
    assert len(output["results"]["vib_freqs"]) == 0
    assert output["results"]["vib_freqs_raw"][0] == 0
    assert output["results"]["vib_freqs_raw"][-1] == 0
    assert output["results"]["vib_freqs"] == []
    assert output["results"]["n_imag"] == 0
    assert output["results"]["imag_vib_freqs"] == []

    assert output["symmetry"]["linear"] is False
    assert output["symmetry"]["rotation_number"] == np.inf
    assert len(output["parameters_thermo"]["vib_freqs"]) == 0
    assert output["results"]["energy"] == -1.0
    assert output["results"]["enthalpy"] == pytest.approx(-0.9357685739989672)
    assert output["results"]["entropy"] == pytest.approx(0.0011292352752446438)
    assert output["results"]["gibbs_energy"] == pytest.approx(-1.2724500713131577)

    atoms = molecule("CH3")
    initial_atoms = deepcopy(atoms)
    output = freq_job(atoms, energy=-10.0, temperature=1000, pressure=20)
    assert output["atoms"] == initial_atoms
    assert len(output["results"]["vib_freqs_raw"]) == 12
    assert len(output["results"]["vib_freqs"]) == 6
    assert output["results"]["vib_energies_raw"][0] == pytest.approx(
        -9.551076713062095e-06
    )
    assert output["results"]["vib_energies_raw"][-1] == pytest.approx(
        0.3880868821616259
    )
    assert output["results"]["vib_energies"][0] == pytest.approx(0.0713506770137291)
    assert output["results"]["vib_energies"][-1] == pytest.approx(0.3880868821616259)
    assert output["results"]["n_imag"] == 0
    assert output["results"]["imag_vib_freqs"] == []

    assert output["parameters_thermo"]["temperature"] == 1000.0
    assert output["parameters_thermo"]["pressure"] == 20.0
    assert output["parameters_thermo"]["sigma"] == 6
    assert output["parameters_thermo"]["spin_multiplicity"] == 2
    assert output["symmetry"]["linear"] is False
    assert output["symmetry"]["rotation_number"] == 6
    assert len(output["parameters_thermo"]["vib_freqs"]) == 6
    assert output["results"]["energy"] == -10.0
    assert output["results"]["enthalpy"] == pytest.approx(-8.749341973959462)
    assert output["results"]["entropy"] == pytest.approx(0.0023506788982171896)
    assert output["results"]["gibbs_energy"] == pytest.approx(-11.100020872176652)
    assert "nid" in output
    assert "dir_name" in output
    assert "nid" in output
    assert "dir_name" in output


def test_unique_workdir(tmpdir):
    from quacc import SETTINGS

    DEFAULT_SETTINGS = SETTINGS.copy()

    SETTINGS.CREATE_UNIQUE_WORKDIR = True
    test_static_job(tmpdir)
    test_relax_job(tmpdir)
    SETTINGS.CREATE_UNIQUE_WORKDIR = DEFAULT_SETTINGS.CREATE_UNIQUE_WORKDIR
