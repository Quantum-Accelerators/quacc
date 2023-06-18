import os
from pathlib import Path
from shutil import rmtree

import pytest
from ase.build import molecule

from quacc.recipes.psi4.core import static_job

try:
    import psi4
except ImportError:
    psi4 = None

FILE_DIR = Path(__file__).resolve().parent


def teardown_module():
    for f in os.listdir(os.getcwd()):
        if f.endswith(".dat"):
            os.remove(f)
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


@pytest.mark.skipif(
    psi4 is None,
    reason="Psi4 must be installed. Try conda install -c psi4 psi4",
)
def test_static_maker():
    atoms = molecule("H2")
    output = static_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["multiplicity"] == 1
    assert output["parameters"]["method"] == "wb97x-v"
    assert output["parameters"]["basis"] == "def2-tzvp"
    assert output["parameters"]["num_threads"] == "max"
    assert output["spin_multiplicity"] == 1
    assert output["charge"] == 0

    output = static_job(
        atoms,
        charge=-2,
        multiplicity=3,
        method="pbe",
        basis="def2-svp",
        swaps={"num_threads": 1, "mem": None, "pop": "regular"},
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["multiplicity"] == 3
    assert output["parameters"]["method"] == "pbe"
    assert output["parameters"]["basis"] == "def2-svp"
    assert output["parameters"]["num_threads"] == 1
    assert output["parameters"]["pop"] == "regular"
    assert "mem" not in output["parameters"]
    assert output["spin_multiplicity"] == 3
    assert output["charge"] == -2
