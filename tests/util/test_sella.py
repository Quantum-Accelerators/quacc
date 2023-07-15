import os
from shutil import rmtree

import pytest
from ase.build import bulk, molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones

from quacc.util.calc import run_ase_opt

try:
    import sella
except ImportError:
    sella = None

CWD = os.getcwd()


def setup_module():
    # Run this test from a fresh directory
    if not os.path.exists("blank_dir"):
        os.mkdir("blank_dir")
    os.chdir("blank_dir")


def teardown_function():
    # Clean up
    os.chdir(CWD)
    for f in os.listdir("."):
        if ".log" in f or ".pckl" in f or ".traj" in f:
            os.remove(f)
    for f in os.listdir(CWD):
        if "quacc-tmp" in f or f == "tmp_dir" or f == "vib" or f == "blank_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


@pytest.mark.skipif(
    sella is None,
    reason="Sella must be installed.",
)
def test_sella():
    from sella.optimize import Sella

    atoms = bulk("Cu") * (2, 1, 1)
    atoms[0].position += 0.1
    atoms.calc = EMT()
    dyn = run_ase_opt(
        atoms,
        optimizer=Sella,
        scratch_dir="test_calc",
        gzip=False,
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.user_internal is False

    atoms = molecule("H2O")
    atoms.calc = LennardJones()
    dyn = run_ase_opt(
        atoms,
        optimizer=Sella,
        scratch_dir="test_calc2",
        gzip=False,
        optimizer_kwargs={"restart": None},
    )
    traj = dyn.traj_atoms
    assert traj[-1].calc.results is not None
    assert dyn.user_internal is True
