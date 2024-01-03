from ase.build import molecule

from quacc.recipes.orca.core import ase_relax_job


def test_ase_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2")

    output = ase_relax_job(atoms, opt_params={"max_steps": 2}, nprocs=1)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["mult"] == 1
    assert (
        output["parameters"]["orcasimpleinput"]
        == "wb97x-d3bj def2-tzvp engrad slowconv normalprint xyzfile"
    )
    assert output.get("trajectory") is not None
    assert len(output["trajectory"]) > 1
    assert output["trajectory"][0] != output["trajectory"][-1]
    assert output.get("attributes") is not None
