from __future__ import annotations

import os
from importlib import util
from pathlib import Path

import pytest
from ase.build import molecule

from quacc.recipes.orca.core import (
    ase_quasi_irc_job,
    ase_relax_job,
    freq_job,
    relax_job,
    static_job,
)

has_fairchem = util.find_spec("fairchem") is not None
has_fairchem_omol = (
    has_fairchem
    and util.find_spec("fairchem.data") is not None
    and util.find_spec("fairchem.data.omol") is not None
)
FILE_DIR = Path(__file__).parent


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = static_job(atoms, charge=0, spin_multiplicity=1, nprocs=1)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-tzvp engrad normalprint wb97m-d4rev xyzfile"
    )
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_static_job_parallel(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = static_job(
        atoms,
        charge=-2,
        spin_multiplicity=3,
        orcasimpleinput=["def2-svp", "#def2-tzvp"],
        orcablocks=["%scf maxiter 300 end"],
        nprocs=2,
    )
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["mult"] == 3
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-svp engrad normalprint wb97m-d4rev xyzfile"
    )
    assert "%scf maxiter 300 end" in output["parameters"]["orcablocks"]


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = relax_job(atoms, charge=0, spin_multiplicity=1, nprocs=2)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-tzvp normalprint opt wb97m-d4rev xyzfile"
    )

    output = relax_job(
        atoms,
        charge=-2,
        spin_multiplicity=3,
        orcasimpleinput=["def2-svp", "#def2-tzvp", "#wb97m-d4rev", "hf"],
        orcablocks=["%scf maxiter 300 end"],
        nprocs=2,
    )
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert (
        output["parameters"]["orcasimpleinput"] == "def2-svp hf normalprint opt xyzfile"
    )
    assert (
        output["parameters"]["orcablocks"] == "%pal nprocs 2 end\n%scf maxiter 300 end"
    )


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_relax_freq_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = relax_job(
        atoms,
        xc="hf",
        basis="def2-svp",
        charge=0,
        spin_multiplicity=1,
        orcasimpleinput=["#normalprint"],
        run_freq=True,
        nprocs=2,
    )
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["orcasimpleinput"] == "def2-svp freq hf opt xyzfile"


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_ase_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = ase_relax_job(atoms, opt_params={"fmax": 0.1}, nprocs=2)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "def2-tzvp engrad normalprint wb97m-d4rev xyzfile"
    )
    assert output["parameters_opt"]["fmax"] == 0.1
    assert output.get("trajectory")
    assert output.get("trajectory_results")


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_ase_relax_job_store(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")
    output = ase_relax_job(
        atoms, opt_params={"store_intermediate_results": True}, nprocs=2
    )
    nsteps = len(output["trajectory"])
    for i in range(nsteps):
        assert f"step{i}" in os.listdir(output["dir_name"])
        assert "orca.xyz.gz" in os.listdir(Path(output["dir_name"], f"step{i}"))


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_freq_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = freq_job(
        atoms,
        xc="hf",
        basis="def2-svp",
        charge=0,
        spin_multiplicity=1,
        orcasimpleinput=["#normalprint"],
        nprocs=2,
    )
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["orcasimpleinput"] == "def2-svp freq hf xyzfile"


@pytest.mark.skipif(os.name == "nt", reason="mpirun not available on Windows")
def test_ase_quasi_irc_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    mode = [[0.0, 0.0, 0.1], [0.0, 0.1, 0.0]]

    output = ase_quasi_irc_job(
        atoms,
        mode,
        charge=0,
        spin_multiplicity=1,
        perturb_magnitude=0.01,
        direction="reverse",
        xc="hf",
        basis="def2-svp",
        orcasimpleinput=["#normalprint"],
        nprocs=2,
    )
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["atoms"] != atoms
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert output["parameters"]["orcasimpleinput"] == "def2-svp engrad hf xyzfile"


@pytest.mark.skipif(not has_fairchem_omol, reason="fairchem not installed")
def test_fairchem_omol(tmp_path, monkeypatch):
    from quacc.recipes.orca.fairchem import omol_static_job

    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")
    output = omol_static_job(atoms)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "ALLPOP def2-TZVPD def2/J DEFGRID3 DIIS EnGrad NONBO NONPA NormalConv NOSOSCF NoUseSym RIJCOSX wB97M-V xyzfile"
    )


@pytest.mark.skipif(not has_fairchem_omol, reason="fairchem not installed")
def test_fairchem_omol_swaps(tmp_path, monkeypatch):
    from quacc.recipes.orca.fairchem import omol_static_job

    monkeypatch.chdir(tmp_path)

    output = omol_static_job(
        molecule("H2"),
        orcasimpleinput=["TightSCF", "#DEFGRID3"],
        orcablocks=["%maxcore 1000"],
        nprocs=2,
    )
    assert (
        output["parameters"]["orcasimpleinput"]
        == "ALLPOP def2-TZVPD def2/J DIIS EnGrad NONBO NONPA NormalConv NOSOSCF NoUseSym RIJCOSX TightSCF wB97M-V xyzfile"
    )
    assert "%maxcore 1000" in output["parameters"]["orcablocks"]
    assert "%pal nprocs 2 end" in output["parameters"]["orcablocks"]
    assert (
        "%scf Convergence Tight maxiter 300 end" in output["parameters"]["orcablocks"]
    )


@pytest.mark.skipif(not has_fairchem_omol, reason="fairchem not installed")
def test_fairchem_omol_open_shell_singlet(tmp_path, monkeypatch):
    from ase import Atoms
    from fairchem.data.omol.orca.calc import Vertical

    from quacc.recipes.orca.fairchem import omol_static_job

    monkeypatch.chdir(tmp_path)

    # Stretched H2: OMol25 ran open-shell singlets with UKS and a 20 degree
    # HOMO-LUMO rotation to break spin symmetry
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 2.5]])
    output = omol_static_job(atoms, vertical=Vertical.Oss, nprocs=2)
    assert "uks" in output["parameters"]["orcasimpleinput"].split()
    assert "%scf rotate {0, 1, 20, 1, 1} end end" in output["parameters"]["orcablocks"]

    # Swaps are still applied alongside the symmetry breaking
    output = omol_static_job(
        atoms, vertical=Vertical.Oss, orcasimpleinput=["TightSCF"], nprocs=2
    )
    assert {"uks", "tightscf"} <= set(output["parameters"]["orcasimpleinput"].split())
    assert "%scf rotate {0, 1, 20, 1, 1} end end" in output["parameters"]["orcablocks"]

    # Triplets are already unrestricted, and the default vertical is unchanged
    output = omol_static_job(
        atoms, vertical=Vertical.Oss, spin_multiplicity=3, nprocs=2
    )
    assert "rotate" not in output["parameters"]["orcablocks"]
    output = omol_static_job(atoms, nprocs=2)
    assert "uks" not in output["parameters"]["orcasimpleinput"].split()
    assert "rotate" not in output["parameters"]["orcablocks"]
