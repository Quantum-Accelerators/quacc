import os

import pytest

from quacc import SETTINGS

DEFAULT_SETTINGS = SETTINGS.copy()


def setup_module():
    import os
    from pathlib import Path

    SETTINGS.WORKFLOW_ENGINE = "local"

    file_dir = Path(__file__).resolve().parent

    with open(file_dir / "mpirun", "w+") as w:
        w.write("")
    os.chmod(file_dir / "mpirun", 0o777)


def teardown_module():
    import os
    from pathlib import Path

    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE
    file_dir = Path(__file__).resolve().parent

    if os.path.exists(file_dir / "mpirun"):
        os.remove(file_dir / "mpirun")


def test_static_job(tmpdir):
    from ase.build import molecule

    from quacc.recipes.orca.core import static_job

    tmpdir.chdir()

    atoms = molecule("H2")

    output = static_job(atoms, 0, 1)
    assert output["natoms"] == len(atoms)
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj def2-tzvp sp slowconv normalprint xyzfile"
    )
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["spin_multiplicity"] == 1
    assert output["charge"] == 0

    output = static_job(
        atoms,
        -2,
        3,
        input_swaps={"def2-svp": True, "def2-tzvp": None},
        block_swaps={"%scf maxiter 300 end": True},
    )
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj sp slowconv normalprint xyzfile def2-svp"
    )
    assert "%scf maxiter 300 end" in output["parameters"]["orcablocks"]


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_relax_job(tmpdir):
    from ase.build import molecule

    from quacc.recipes.orca.core import relax_job

    tmpdir.chdir()

    atoms = molecule("H2")

    output = relax_job(atoms, 0, 1)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj def2-tzvp opt slowconv normalprint xyzfile"
    )
    assert output["trajectory"][0] != output["trajectory"][-1]

    output = relax_job(
        atoms,
        -2,
        3,
        input_swaps={
            "hf": True,
            "wb97x-d3bj": None,
            "def2-svp": True,
            "def2-tzvp": None,
        },
        block_swaps={"%scf maxiter 300 end": True},
    )
    assert output["natoms"] == len(atoms)
    assert (
        output["parameters"]["orcasimpleinput"]
        == "opt slowconv normalprint xyzfile hf def2-svp"
    )
    assert "%scf maxiter 300 end" in output["parameters"]["orcablocks"]
    assert "trajectory" in output
    assert len(output["trajectory"]) > 1


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_mpi_run(tmpdir, monkeypatch):
    import multiprocessing
    from pathlib import Path

    from ase.build import molecule

    from quacc.recipes.orca.core import relax_job, static_job

    file_dir = Path(__file__).resolve().parent
    tmpdir.chdir()
    monkeypatch.setenv("PATH", file_dir)

    atoms = molecule("H2")
    output = static_job(atoms, 0, 1)
    nprocs = multiprocessing.cpu_count()
    assert f"%pal nprocs {nprocs} end" in output["parameters"]["orcablocks"]

    output = relax_job(atoms, 0, 1)
    nprocs = multiprocessing.cpu_count()
    assert f"%pal nprocs {nprocs} end" in output["parameters"]["orcablocks"]
