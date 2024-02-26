from ase.build import bulk

from quacc.recipes.emt.elastic import bulk_to_deformations_flow


def test_elastic_jobs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    outputs = bulk_to_deformations_flow(atoms, run_static=False)
    assert [output["parameters"]["asap_cutoff"] is False for output in outputs]
    assert [output["name"] == "EMT Relax" for output in outputs]
