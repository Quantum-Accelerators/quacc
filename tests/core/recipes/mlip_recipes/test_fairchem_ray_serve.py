from __future__ import annotations

from importlib.util import find_spec

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("ray")
pytest.importorskip("fairchem")
pytest.importorskip("fairchem.core")

if find_spec("fairchem"):
    from huggingface_hub.utils._auth import get_token

    if not get_token():
        pytest.skip(
            "HF_TOKEN required for fairchem ray-serve tests", allow_module_level=True
        )

import numpy as np
from ase.build import bulk

from quacc import get_settings
from quacc.recipes.mlip._base import pick_calculator
from quacc.recipes.mlip.core import relax_job, static_job
from quacc.recipes.mlip.elastic import elastic_tensor_flow


@pytest.fixture(scope="module")
def ray_serve_cluster():
    """Spin up a local Ray cluster + multiplexed fairchem serve deployment.

    The deployment is named ``multiplexed-predict-server`` to match what
    ``quacc.recipes.mlip._base`` connects to when
    ``FAIRCHEM_RAY_SERVE_BATCHING`` is enabled.
    """
    import ray
    from fairchem.core.components.batch_server import (
        setup_multiplexed_batch_predict_server,
    )
    from ray import serve

    started_ray = False
    if not ray.is_initialized():
        ray.init(
            include_dashboard=False,
            ignore_reinit_error=True,
            log_to_driver=False,
            configure_logging=False,
        )
        started_ray = True

    handle = setup_multiplexed_batch_predict_server(
        deployment_config={"num_replicas": 1, "ray_actor_options": {"num_cpus": 1}},
        deployment_name="multiplexed-predict-server",
        route_prefix="/multiplex-predict",
    )
    try:
        yield handle
    finally:
        try:
            serve.shutdown()
        finally:
            if started_ray and ray.is_initialized():
                ray.shutdown()


@pytest.fixture
def _batching(monkeypatch, request, ray_serve_cluster):
    settings = get_settings()
    monkeypatch.setattr(
        settings, "FAIRCHEM_RAY_SERVE_BATCHING", request.param, raising=False
    )
    # pick_calculator is wrapped by freezeargs(lru_cache(...)); clear the
    # underlying lru_cache so each parametrize value rebuilds the calc
    # under the current FAIRCHEM_RAY_SERVE_BATCHING setting.
    pick_calculator.__wrapped__.cache_clear()
    yield request.param
    pick_calculator.__wrapped__.cache_clear()


def _assert_calc_path(batching):
    """Verify the calculator built during the test actually used the
    expected code path (Ray Serve batching vs. local inference)."""
    from fairchem.core.units.mlip_unit.predict import BatchServerPredictUnit

    calc = pick_calculator(
        library="fairchem", name_or_path="uma-s-1p1", task_name="omat"
    )
    used_batching = isinstance(getattr(calc, "predictor", None), BatchServerPredictUnit)
    assert used_batching is batching, (
        f"Expected FAIRCHEM_RAY_SERVE_BATCHING={batching} path, "
        f"but BatchServerPredictUnit usage was {used_batching}."
    )


@pytest.mark.parametrize("_batching", [False, True], indirect=True)
def test_static_job(tmp_path, monkeypatch, _batching):
    monkeypatch.chdir(tmp_path)
    torch.set_default_dtype(torch.float32)

    atoms = bulk("Cu")
    output = static_job(
        atoms, library="fairchem", name_or_path="uma-s-1p1", task_name="omat"
    )
    assert output["results"]["energy"] == pytest.approx(-3.7501682869643735, rel=1e-3)
    assert np.shape(output["results"]["forces"]) == (1, 3)
    assert output["atoms"] == atoms
    _assert_calc_path(_batching)


@pytest.mark.parametrize("_batching", [False, True], indirect=True)
def test_relax_job(tmp_path, monkeypatch, _batching):
    monkeypatch.chdir(tmp_path)
    torch.set_default_dtype(torch.float32)

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(
        atoms, library="fairchem", name_or_path="uma-s-1p1", task_name="omat"
    )
    assert output["results"]["energy"] == pytest.approx(-30.001143639922756, rel=1e-3)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() == pytest.approx(atoms.get_volume())
    _assert_calc_path(_batching)


@pytest.mark.parametrize("_batching", [False, True], indirect=True)
def test_relax_cell_job(tmp_path, monkeypatch, _batching):
    monkeypatch.chdir(tmp_path)
    torch.set_default_dtype(torch.float32)

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].position += 0.1
    output = relax_job(
        atoms,
        library="fairchem",
        relax_cell=True,
        name_or_path="uma-s-1p1",
        task_name="omat",
    )
    assert output["results"]["energy"] == pytest.approx(-30.005004590392726, rel=1e-3)
    assert np.shape(output["results"]["forces"]) == (8, 3)
    assert output["atoms"] != atoms
    assert output["atoms"].get_volume() != pytest.approx(atoms.get_volume())
    _assert_calc_path(_batching)


@pytest.mark.parametrize("_batching", [False, True], indirect=True)
def test_elastic_jobs(tmp_path, monkeypatch, _batching):
    monkeypatch.chdir(tmp_path)
    torch.set_default_dtype(torch.float32)

    atoms = bulk("Cu")
    outputs = elastic_tensor_flow(
        atoms,
        run_static=False,
        pre_relax=True,
        job_params={
            "all": {
                "library": "fairchem",
                "name_or_path": "uma-s-1p1",
                "task_name": "omat",
            },
            "relax_job": {"opt_params": {"fmax": 0.01}},
        },
    )
    assert outputs["deformed_results"][0]["atoms"].get_volume() != pytest.approx(
        atoms.get_volume()
    )
    assert outputs["undeformed_result"]["results"]["stress"] == pytest.approx(
        0, abs=1e-2
    )
    assert outputs["elasticity_doc"].bulk_modulus.voigt == pytest.approx(151.367, abs=2)
    _assert_calc_path(_batching)
