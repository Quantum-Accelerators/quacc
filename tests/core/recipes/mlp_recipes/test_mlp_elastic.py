from __future__ import annotations

import pytest
from ase.build import bulk

from quacc.recipes.mlp.elastic import bulk_to_deformations_flow


def test_elastic_jobs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")

    outputs = bulk_to_deformations_flow(
        atoms,
        run_static=False,
        pre_relax=True,
        job_params={"all": {"method": "sevennet"}},
    )
    assert outputs["deformed_results"][0]["atoms"].get_volume() != pytest.approx(
        atoms.get_volume()
    )
    assert outputs["undeformed_result"]["results"]["stress"] == pytest.approx(
        0, abs=1e-2
    )
    assert outputs["elasticity_doc"].bulk_modulus.voigt == pytest.approx(
        143.771, abs=1e-2
    )
    for output in outputs["deformed_results"]:
        assert output["nelements"] == 1
        assert output["nsites"] == 1
    assert len(outputs["deformed_results"]) == 24

    outputs = bulk_to_deformations_flow(
        atoms,
        run_static=True,
        pre_relax=True,
        job_params={"all": {"method": "sevennet"}},
    )
    assert outputs["deformed_results"][0]["atoms"].get_volume() != pytest.approx(
        atoms.get_volume()
    )

    for output in outputs["deformed_results"]:
        assert output["nelements"] == 1
        assert output["nsites"] == 1
    assert len(outputs["deformed_results"]) == 24
