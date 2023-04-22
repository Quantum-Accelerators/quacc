import os
from shutil import rmtree

import numpy as np
import pytest
from ase.build import bulk, molecule

try:
    from xtb.ase.calculator import XTB
except ImportError:
    XTB = None
from quacc.recipes.xtb.core import relax_job, static_job, thermo_job


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
    XTB is None,
    reason="xTB-python must be installed. Try conda install -c conda-forge xtb-python",
)
def test_static_Job():
    atoms = molecule("H2O")
    output = static_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.96777587302995)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    output = static_job(atoms, method="GFN1-xTB")
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["results"]["energy"] == pytest.approx(-156.96750578831137)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    output = static_job(atoms, method="GFN-FF")
    assert output["parameters"]["method"] == "GFN-FF"
    assert output["results"]["energy"] == pytest.approx(-8.912667188932252)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    atoms = bulk("Cu")
    output = static_job(atoms, method="GFN1-xTB")
    assert output["results"]["energy"] == pytest.approx(-119.77643232313169)
    assert np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


@pytest.mark.skipif(
    XTB is None,
    reason="xTB-python must be installed. Try conda install -c conda-forge xtb-python",
)
def test_relax_Job():
    atoms = molecule("H2O")
    output = relax_job(atoms)
    assert output["spin_multiplicity"] == 1
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["method"] == "GFN2-xTB"
    assert output["results"]["energy"] == pytest.approx(-137.9764670127011)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    output = relax_job(atoms, method="GFN1-xTB")
    assert output["parameters"]["method"] == "GFN1-xTB"
    assert output["results"]["energy"] == pytest.approx(-156.9763496338962)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    output = relax_job(atoms, method="GFN-FF")
    assert output["parameters"]["method"] == "GFN-FF"
    assert output["results"]["energy"] == pytest.approx(-8.915974748299963)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())

    atoms = molecule("H2O")
    atoms.center(vacuum=5)
    atoms.pbc = True
    atoms[0].position += 0.01
    output = relax_job(atoms, method="GFN1-xTB")
    assert output["results"]["energy"] == pytest.approx(-156.97441169886613)
    assert not np.array_equal(output["atoms"].get_positions(), atoms.get_positions())


@pytest.mark.skipif(
    XTB is None,
    reason="xTB-python must be installed. Try conda install -c conda-forge xtb-python",
)
def test_thermo_job():
    atoms = molecule("H2O")
    output = thermo_job(atoms)
    assert output["atoms"] == atoms
    assert output["results"]["n_imag"] == 0
    assert len(output["results"]["frequencies"]) == 9
    assert len(output["results"]["true_frequencies"]) == 3
    assert output["results"]["true_frequencies"][-1] == pytest.approx(3526.945468014458)
    assert output["results"]["geometry"] == "nonlinear"
    assert output["results"]["energy"] == 0.0
    assert output["results"]["enthalpy"] == pytest.approx(0.637581401404518)
    assert output["results"]["entropy"] == pytest.approx(0.0019584993671715764)
    assert output["results"]["gibbs_energy"] == pytest.approx(0.05365481508231251)

    atoms = molecule("O2")
    output = thermo_job(atoms, energy=-100.0, temperature=200, pressure=2.0)
    assert output["atoms"] == atoms
    assert output["results"]["n_imag"] == 0
    assert len(output["results"]["true_frequencies"]) == 1
    assert output["results"]["true_frequencies"][-1] == pytest.approx(1449.82397338371)
    assert output["results"]["geometry"] == "linear"
    assert output["results"]["pointgroup"] == "D*h"
    assert output["results"]["energy"] == -100.0
    assert output["results"]["enthalpy"] == pytest.approx(-99.84979574721257)
    assert output["results"]["entropy"] == pytest.approx(0.001915519747865423)
    assert output["results"]["gibbs_energy"] == pytest.approx(-100.23289969678565)

    atoms = molecule("H")
    output = thermo_job(atoms, energy=-1.0)
    assert output["atoms"] == atoms
    assert output["results"]["n_imag"] == 0
    assert output["results"]["geometry"] == "monatomic"
    assert len(output["results"]["true_frequencies"]) == 0
    assert output["results"]["energy"] == -1.0
    assert output["results"]["enthalpy"] == pytest.approx(-0.9357685739989672)
