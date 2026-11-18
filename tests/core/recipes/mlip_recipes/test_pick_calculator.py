"""Tests for the caching of ``pick_calculator``."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace

import pytest

from quacc.recipes.mlip._base import pick_calculator


@pytest.fixture
def fake_matcalc(monkeypatch):
    calls = []

    def load_fp(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(parameters={})

    module = ModuleType("matcalc")
    module.__version__ = "0.0.0"
    module.load_fp = load_fp
    monkeypatch.setitem(sys.modules, "matcalc", module)

    pick_calculator.__wrapped__.cache_clear()
    yield calls
    pick_calculator.__wrapped__.cache_clear()


def test_pick_calculator_nested_kwargs(fake_matcalc):
    overrides = {"backbone": {"always_use_pbc": False, "cutoffs": [5.0, 6.0]}}
    calc = pick_calculator("matcalc", name="model", overrides=overrides, workers=[1])

    # The library receives the kwargs exactly as they were passed
    assert fake_matcalc == [{"name": "model", "overrides": overrides, "workers": [1]}]
    assert type(fake_matcalc[0]["overrides"]) is dict
    assert type(fake_matcalc[0]["overrides"]["backbone"]["cutoffs"]) is list

    # Identical (but distinct) kwargs reuse the cached calculator
    same_calc = pick_calculator(
        "matcalc",
        name="model",
        overrides={"backbone": {"always_use_pbc": False, "cutoffs": [5.0, 6.0]}},
        workers=[1],
    )
    assert same_calc is calc
    assert len(fake_matcalc) == 1

    # Different nested kwargs build a new calculator
    pick_calculator(
        "matcalc", name="model", overrides={"backbone": {"always_use_pbc": True}}
    )
    assert len(fake_matcalc) == 2


def test_freeze_thaw_roundtrip():
    from quacc.recipes.mlip._base import freeze, thaw

    value = {"a": [1, {"b": (2, [3])}], "c": "d"}
    frozen = freeze(value)
    hash(frozen)
    assert thaw(frozen) == value
    assert type(thaw(frozen)["a"]) is list
    assert type(thaw(frozen)["a"][1]["b"]) is tuple
