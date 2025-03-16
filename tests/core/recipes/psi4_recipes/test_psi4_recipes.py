from __future__ import annotations

import pytest

pytest.importorskip("psi4")

from ase.build import molecule

from quacc.recipes.psi4.core import static_job


def test_static(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")
    output = static_job(atoms, charge=0, spin_multiplicity=1)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["multiplicity"] == 1
    assert output["parameters"]["method"] == "wb97x-v"
    assert output["parameters"]["basis"] == "def2-tzvp"
    assert output["parameters"]["num_threads"] == "max"

    output = static_job(
        atoms,
        charge=-2,
        spin_multiplicity=3,
        method="pbe",
        basis="def2-svp",
        num_threads=1,
        mem=None,
        pop="regular",
    )
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == -2
    assert output["parameters"]["multiplicity"] == 3
    assert output["parameters"]["method"] == "pbe"
    assert output["parameters"]["basis"] == "def2-svp"
    assert output["parameters"]["pop"] == "regular"
    assert output["parameters"]["num_threads"] == 1
    assert "mem" not in output["parameters"]
