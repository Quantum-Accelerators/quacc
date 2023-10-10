import pytest
from ase.build import molecule

from quacc import SETTINGS

psi4 = pytest.importorskip("psi4")

pytestmark = pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE != "local",
    reason="Need to use local as workflow manager to run this test.",
)


def test_static(tmpdir):
    from quacc.recipes.psi4.core import static_job

    tmpdir.chdir()

    atoms = molecule("H2")
    output = static_job(atoms, 0, 1)
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
        -2,
        3,
        method="pbe",
        basis="def2-svp",
        calc_swaps={"num_threads": 1, "mem": None, "pop": "regular"},
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
