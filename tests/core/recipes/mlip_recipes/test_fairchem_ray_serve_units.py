"""Unit tests for the Ray-Serve branch of ``pick_calculator``.

These tests stub out the heavy fairchem and Ray Serve calls so they can
run without HF_TOKEN, without a real Ray cluster, and without actually
loading any UMA checkpoint. They cover the small but easy-to-miss
branches in ``pick_calculator``: the two fallback paths (Ray missing,
Ray uninitialized) and the alternative model-identifier kwarg paths
(``model_id``/``checkpoint``).
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

import pytest

pytest.importorskip("fairchem")
pytest.importorskip("fairchem.core")
pytest.importorskip("ray")

from fairchem.core.units.mlip_unit.api.model_spec import ModelSpec

from quacc import get_settings
from quacc.recipes.mlip._base import pick_calculator


@pytest.fixture
def enable_batching(monkeypatch):
    settings = get_settings()
    monkeypatch.setattr(settings, "FAIRCHEM_RAY_SERVE_BATCHING", True, raising=False)
    pick_calculator.__wrapped__.cache_clear()
    yield
    pick_calculator.__wrapped__.cache_clear()


def _stub_calc(**_kwargs):
    """Return a sentinel object so we can assert the local-fallback path
    constructed a calculator without actually downloading a checkpoint."""
    return SimpleNamespace(parameters={})


@pytest.mark.usefixtures("enable_batching")
def test_falls_back_when_ray_not_initialized(monkeypatch, caplog):
    import ray

    monkeypatch.setattr(ray, "is_initialized", lambda: False)
    with (
        patch(
            "fairchem.core.FAIRChemCalculator.from_model_checkpoint",
            side_effect=_stub_calc,
        ) as mock_local,
        caplog.at_level("WARNING"),
    ):
        pick_calculator(library="fairchem", name_or_path="uma-s-1p1")
    mock_local.assert_called_once()
    assert "Ray is not initialized" in caplog.text


@pytest.mark.usefixtures("enable_batching")
def test_falls_back_when_ray_not_installed(monkeypatch, caplog):
    import builtins

    real_import = builtins.__import__

    def _no_ray(name, *args, **kwargs):
        if name == "ray":
            raise ImportError("ray not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _no_ray)
    with (
        patch(
            "fairchem.core.FAIRChemCalculator.from_model_checkpoint",
            side_effect=_stub_calc,
        ) as mock_local,
        caplog.at_level("WARNING"),
    ):
        pick_calculator(library="fairchem", name_or_path="uma-s-1p1-fallback")
    mock_local.assert_called_once()
    assert "Ray is not installed" in caplog.text


def _stub_predict_unit(**_kwargs):
    return _kwargs


@pytest.fixture
def stub_serve_unit():
    """Replace ``BatchServerPredictUnit.from_deployment_connection_info``
    so we don't need a live Ray Serve deployment."""
    captured = {}

    def _capture(**kwargs):
        captured.update(kwargs)
        return SimpleNamespace(
            dataset_to_tasks={"omat": [SimpleNamespace(property="energy")]},
            inference_settings=SimpleNamespace(
                external_graph_gen=False, base_precision_dtype="float32"
            ),
        )

    with (
        patch(
            "fairchem.core.units.mlip_unit.predict.BatchServerPredictUnit.from_deployment_connection_info",
            side_effect=_capture,
        ),
        patch(
            "fairchem.core.FAIRChemCalculator",
            side_effect=lambda **k: SimpleNamespace(
                parameters={}, predictor=k.get("predict_unit")
            ),
        ),
    ):
        yield captured


@pytest.mark.usefixtures("enable_batching")
def test_serve_branch_uses_name_or_path(monkeypatch, stub_serve_unit):
    import ray

    monkeypatch.setattr(ray, "is_initialized", lambda: True)
    pick_calculator(
        library="fairchem", name_or_path="my/local/ckpt.pt", task_name="oc20"
    )
    assert stub_serve_unit["model_spec"].checkpoint == "my/local/ckpt.pt"
    assert stub_serve_unit["deployment_name"] == "multiplexed-predict-server"


@pytest.mark.usefixtures("enable_batching")
def test_serve_branch_uses_model_id(monkeypatch, stub_serve_unit):
    import ray

    monkeypatch.setattr(ray, "is_initialized", lambda: True)
    pick_calculator(
        library="fairchem",
        model_id="uma-s-2",
        inference_settings="default",
        task_name="omat",
    )
    assert stub_serve_unit["model_spec"].checkpoint == "uma-s-2"


@pytest.mark.usefixtures("enable_batching")
def test_serve_branch_default_checkpoint(monkeypatch, stub_serve_unit):
    import ray

    monkeypatch.setattr(ray, "is_initialized", lambda: True)
    # Neither name_or_path, model_id, nor checkpoint provided → default
    pick_calculator(library="fairchem", task_name="omat")
    assert stub_serve_unit["model_spec"].checkpoint == "uma-s-1p1"


@pytest.mark.usefixtures("enable_batching")
def test_serve_branch_builds_model_spec(monkeypatch, stub_serve_unit):
    """Model-loading kwargs are represented by FairChem's ``ModelSpec``."""
    import ray

    monkeypatch.setattr(ray, "is_initialized", lambda: True)
    pick_calculator(
        library="fairchem",
        name_or_path="uma-s-1p1",
        task_name="omat",
        device="cpu",
        overrides={"foo": 1},
        seed=42,
    )
    model_spec = stub_serve_unit["model_spec"]
    assert isinstance(model_spec, ModelSpec)
    assert model_spec.checkpoint == "uma-s-1p1"
    assert model_spec.device == "cpu"
    assert model_spec.overrides == {"foo": 1}
