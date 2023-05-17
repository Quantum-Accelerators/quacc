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
    assert output["atoms"] == atoms
    assert output["results"]["n_imag"] == 0
    assert len(output["results"]["frequencies"]) == 9
    assert len(output["results"]["true_frequencies"]) == 3
    assert output["results"]["true_frequencies"][-1] == pytest.approx(
        3526.9940431751647
    )
    assert output["results"]["geometry"] == "nonlinear"
    assert output["results"]["energy"] == 0.0
    assert output["results"]["enthalpy"] == pytest.approx(0.6375973622705744)
    assert output["results"]["entropy"] == pytest.approx(0.0019584992229988523)
    assert output["results"]["gibbs_energy"] == pytest.approx(0.05367081893346748)

    atoms = molecule("O2")
    output = thermo_job(atoms, energy=-100.0, temperature=200, pressure=2.0)
    assert output["atoms"] == atoms
    assert output["results"]["n_imag"] == 0
    assert len(output["results"]["true_frequencies"]) == 1
    assert output["results"]["true_frequencies"][-1] == pytest.approx(
        1449.8292119293476
    )
    assert output["results"]["geometry"] == "linear"
    assert output["results"]["pointgroup"] == "D*h"
    assert output["results"]["energy"] == -100.0
    assert output["results"]["enthalpy"] == pytest.approx(-99.84979543613858)
    assert output["results"]["entropy"] == pytest.approx(0.0019155197468650228)
    assert output["results"]["gibbs_energy"] == pytest.approx(-100.23289938549244)

    atoms = molecule("H")
    output = thermo_job(atoms, energy=-1.0)
    assert output["atoms"] == atoms
    assert output["results"]["n_imag"] == 0
    assert output["results"]["geometry"] == "monatomic"
    assert len(output["results"]["true_frequencies"]) == 0
    assert output["results"]["energy"] == -1.0
    assert output["results"]["enthalpy"] == pytest.approx(-0.9357685739989672)
