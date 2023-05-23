import os
from shutil import rmtree

import numpy as np
import pytest
from ase.build import molecule

try:
    from tblite.ase import TBLite
except ImportError:
    TBLite = None
from quacc.recipes.tblite.core import relax_job, static_job, thermo_job


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
    reason="tblite must be installed. Try pip install tblite[ase]",
)
def test_static_Job():
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
    reason="tblite must be installed. Try pip install tblite[ase]",
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
    reason="tblite must be installed. Try pip install tblite[ase]",
)
def test_thermo_job():
    atoms = molecule("H2O")
    output = thermo_job(atoms)
    assert output["vib"]["atoms"] == atoms
    assert len(output["vib"]["results"]["vib_freqs"]) == 9
    assert len(output["vib"]["results"]["true_vib_freqs"]) == 3
    assert output["vib"]["results"]["vib_freqs"][0] == pytest.approx(
        -0.10864429415434408
    )
    assert output["vib"]["results"]["vib_freqs"][-1] == pytest.approx(
        3526.9940431752034
    )
    assert output["vib"]["results"]["true_vib_freqs"][0] == pytest.approx(
        1586.623114694335
    )
    assert output["vib"]["results"]["true_vib_freqs"][-1] == pytest.approx(
        3526.9940431752034
    )
    assert output["vib"]["results"]["n_imag"] == 0
    assert output["vib"]["results"]["imag_vib_freqs"] == []

    assert output["thermo"]["atoms"] == atoms
    assert output["thermo"]["symmetry"]["point_group"] == "C2v"
    assert output["thermo"]["symmetry"]["rotation_number"] == 2
    assert output["thermo"]["symmetry"]["linear"] is False
    assert len(output["thermo"]["results"]["vib_freqs"]) == 3
    assert output["vib"]["results"]["true_vib_freqs"][0] == pytest.approx(
        1586.623114694335
    )
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

    # atoms = molecule("O2")
    # output = thermo_job(atoms, energy=-100.0, temperature=200, pressure=2.0)
    # assert output["vib"]["atoms"] == atoms
    # assert output["vib"]["results"]["n_imag"] == 0

    # assert output["thermo"]["atoms"] == atoms
    # assert len(output["thermo"]["results"]["vib_freqs"]) == 1
    # assert output["thermo"]["results"]["vib_freqs"][-1] == pytest.approx(
    #     1449.8291293360887
    # )
    # assert output["thermo"]["symmetry"]["linear"] is True
    # assert output["thermo"]["symmetry"]["point_group"] == "D*h"
    # assert output["thermo"]["symmetry"]["rotation_number"] == 2
    # assert output["thermo"]["results"]["energy"] == -100.0
    # assert output["thermo"]["results"]["enthalpy"] == pytest.approx(-99.849795427762328)
    # assert output["thermo"]["results"]["entropy"] == pytest.approx(
    #     0.0019504600146485497
    # )
    # assert output["thermo"]["results"]["gibbs_energy"] == pytest.approx(
    #     -100.23988743069202
    # )

    # atoms = molecule("H")
    # output = thermo_job(atoms, energy=-1.0)
    # assert output["vib"]["atoms"] == atoms
    # assert output["vib"]["results"]["n_imag"] == 0

    # assert output["thermo"]["atoms"] == atoms
    # assert output["thermo"]["symmetry"]["linear"] is False
    # assert output["thermo"]["symmetry"]["rotation_number"] == np.inf
    # assert len(output["thermo"]["results"]["vib_freqs"]) == 0
    # assert output["thermo"]["results"]["energy"] == -1.0
    # assert output["thermo"]["results"]["enthalpy"] == pytest.approx(-0.9357685739989672)
