from __future__ import annotations

import pytest
from ase.build import bulk

from quacc.recipes.emt.elastic import bulk_to_elastic_tensor_flow


def test_elastic_jobs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    outputs = bulk_to_elastic_tensor_flow(atoms, run_static=False)
    assert outputs["deformed_results"][0]["atoms"].get_volume() != pytest.approx(
        atoms.get_volume()
    )
    assert outputs["elasticity_doc"].bulk_modulus.voigt == pytest.approx(134.579)
    for output in outputs["deformed_results"]:
        assert output["parameters"]["asap_cutoff"] is False
        assert output["name"] == "EMT Relax"
        assert output["nelements"] == 1
        assert output["nsites"] == 1
    assert len(outputs["deformed_results"]) == 24

    outputs = bulk_to_elastic_tensor_flow(
        atoms, run_static=True, job_params={"static_job": {"asap_cutoff": True}}
    )
    assert outputs["deformed_results"][0]["atoms"].get_volume() != pytest.approx(
        atoms.get_volume()
    )
    assert outputs["deformed_results"][0]["atoms"].get_volume() != pytest.approx(
        atoms.get_volume()
    )
    for output in outputs["deformed_results"]:
        assert output["parameters"]["asap_cutoff"] is True
        assert output["name"] == "EMT Static"
        assert output["nelements"] == 1
        assert output["nsites"] == 1
    assert len(outputs["deformed_results"]) == 24
