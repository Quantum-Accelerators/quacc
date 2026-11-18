from __future__ import annotations

from types import SimpleNamespace

import pytest

pytest.importorskip("torch_sim")

from quacc.recipes.torchsim._base import process_trajectory_reporter_dict
from quacc.schemas.torchsim import PROPERTY_FN_REGISTRY, PropertyFn


def test_property_fns_are_registered():
    for prop in PropertyFn:
        assert prop in PROPERTY_FN_REGISTRY
    assert PropertyFn.ENERGY == "potential_energy"


def test_trajectory_reporter_accepts_property_fns(tmp_path):
    runner = SimpleNamespace(tmpdir=tmp_path)
    _, details = process_trajectory_reporter_dict(
        {"prop_calculators": {1: [PropertyFn.ENERGY, PropertyFn.TEMPERATURE]}},
        runner,  # type: ignore[arg-type]
        n_systems=1,
    )
    assert details["prop_calculators"] == {1: ["potential_energy", "temperature"]}
