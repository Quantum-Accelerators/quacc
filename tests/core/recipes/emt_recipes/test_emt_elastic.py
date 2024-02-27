from ase.build import bulk
from numpy.testing import assert_equal

from quacc.recipes.emt.elastic import bulk_to_deformations_flow


def test_elastic_jobs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    outputs = bulk_to_deformations_flow(atoms, run_static=True)
    for output in outputs:
        assert output["parameters"]["asap_cutoff"] is False
        assert output["name"] == "EMT Relax"
        assert output["nelements"] == 1
        assert output["nsites"] == 1
        assert len(outputs) == 24
