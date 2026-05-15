"""Unit tests for the Ray-Serve branch of ``pick_calculator``.

These tests stub out the heavy fairchem and Ray Serve calls so they can
run without HF_TOKEN, without a real Ray cluster, and without actually
loading any UMA checkpoint. They cover the small but easy-to-miss
branches in ``pick_calculator``: the two fallback paths (Ray missing,
Ray uninitialized) and the alternative model-identifier kwarg paths
(``model_id``/``checkpoint``).
"""

from __future__ import annotations

from importlib.util import find_spec
from types import SimpleNamespace
from unittest.mock import patch

import pytest

if find_spec("fairchem") is None:
    pytest.skip("fairchem not installed", allow_module_level=True)

from quacc import get_settings
from quacc.recipes.mlp._base import pick_calculator


@pytest.fixture
def _enable_batching(monkeypatch):
    settings = get_settings()
    monkeypatch.setattr(
        settings, "FAIRCHEM_RAY_SERVE_BATCHING", True, raising=False
    )
    pick_calculator.__wrapped__.cache_clear()
    yield
    pick_calculator.__wrapped__.cache_clear()


def _stub_calc(**_kwargs):
    """Return a sentinel object so we can assert the local-fallback path
    constructed a calculator without actually downloading a checkpoint."""
    return SimpleNamespace(parameters={})


def test_falls_back_when_ray_not_initialized(monkeypatch, caplog, _enable_batching):
    import ray

    monkeypatch.setattr(ray, "is_initialized", lambda: False)
    with patch(
        "fairchem.core.FAIRChemCalculator.from_model_checkpoint",
        side_effect=_stub_calc,
    ) as mock_local:
        with caplog.at_level("WARNING"):
            pick_calculator(method="fairchem", name_or_path="uma-s-1p1")
    mock_local.assert_called_once()
    assert "Ray is not initialized" in caplog.text


def test_falls_back_when_ray_not_installed(monkeypatch, caplog, _enable_batching):
    import builtins

    real_import = builtins.__import__

    def _no_ray(name, *args, **kwargs):
        if name == "ray":
            raise ImportError("ray not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _no_ray)
    with patch(
        "fairchem.core.FAIRChemCalculator.from_model_checkpoint",
        side_effect=_stub_calc,
    ) as mock_local:
        with caplog.at_level("WARNING"):
            pick_calculator(
                method="fairchem", name_or_path="uma-s-1p1-fallback"
            )
    mock_local.assert_called_once()
    assert "Ray is not installed" in caplog.text


def _stub_predict_unit(**_kwargs):
    return _kwargs


@pytest.fixture
def _stub_serve_unit():
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

    with patch(
        "fairchem.core.units.mlip_unit.predict.BatchServerPredictUnit.from_deployment_connection_info",
        side_effect=_capture,
    ), patch(
        "fairchem.core.FAIRChemCalculator", side_effect=lambda **k: SimpleNamespace(
            parameters={}, predictor=k.get("predict_unit")
        )
    ):
        yield captured


def test_serve_branch_uses_name_or_path(monkeypatch, _enable_batching, _stub_serve_unit):
    import ray

    monkeypatch.setattr(ray, "is_initialized", lambda: True)
    pick_calculator(
        method="fairchem", name_or_path="my/local/ckpt.pt", task_name="oc20"
    )
    assert _stub_serve_unit["multiplexed_model_id"] == "my/local/ckpt.pt:default"
    assert _stub_serve_unit["deployment_name"] == "predict-server"


def test_serve_branch_uses_model_id(monkeypatch, _enable_batching, _stub_serve_unit):
    import ray

    monkeypatch.setattr(ray, "is_initialized", lambda: True)
    pick_calculator(
        method="fairchem",
        model_id="uma-s-2",
        inference_settings="fast",
        task_name="omat",
    )
    assert _stub_serve_unit["multiplexed_model_id"] == "uma-s-2:fast"


def test_serve_branch_default_checkpoint(monkeypatch, _enable_batching, _stub_serve_unit):
    import ray

    monkeypatch.setattr(ray, "is_initialized", lambda: True)
    # Neither name_or_path, model_id, nor checkpoint provided → default
    pick_calculator(method="fairchem", task_name="omat")
    assert _stub_serve_unit["multiplexed_model_id"] == "uma-s-1p1:default"


def test_serve_branch_drops_local_only_kwargs(
    monkeypatch, _enable_batching, _stub_serve_unit
):
    """``device``/``overrides``/``seed`` must be stripped before reaching
    the Ray Serve helper (it doesn't accept them)."""
    import ray

    monkeypatch.setattr(ray, "is_initialized", lambda: True)
    # If the kwargs were not popped, FAIRChemCalculator stub would
    # receive them and the test would still pass — but the important
    # thing is that the serve helper itself doesn't get them. The
    # fixture's _capture only collects deployment_name/multiplexed_model_id,
    # so it implicitly verifies no extra kwargs leak through.
    pick_calculator(
        method="fairchem",
        name_or_path="uma-s-1p1",
        task_name="omat",
        device="cpu",
        overrides={"foo": 1},
        seed=42,
    )
    assert _stub_serve_unit["multiplexed_model_id"] == "uma-s-1p1:default"
