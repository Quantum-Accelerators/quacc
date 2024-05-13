from __future__ import annotations

import pytest
from ase.build import bulk

from quacc.recipes.emt.elastic import bulk_to_deformations_flow


def test_elastic_jobs():
    atoms = bulk("Cu")

    outputs = bulk_to_deformations_flow(atoms, run_static=False)
    assert outputs[0]["atoms"].get_volume() != pytest.approx(atoms.get_volume())
    for output in outputs:
        assert output["parameters"]["asap_cutoff"] is False
        assert output["name"] == "EMT Relax"
        assert output["nelements"] == 1
        assert output["nsites"] == 1
        assert len(outputs) == 24

    outputs = bulk_to_deformations_flow(
        atoms, run_static=True, job_params={"static_job": {"asap_cutoff": True}}
    )
    assert outputs[0]["atoms"].get_volume() != pytest.approx(atoms.get_volume())
    for output in outputs:
        assert output["parameters"]["asap_cutoff"] is True
        assert output["name"] == "EMT Static"
        assert output["nelements"] == 1
        assert output["nsites"] == 1
        assert len(outputs) == 24
