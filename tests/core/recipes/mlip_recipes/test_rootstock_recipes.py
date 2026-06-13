"""Mocked unit tests for rootstock MLIP integration in quacc."""

from __future__ import annotations

import sys
import types

import numpy as np
import pytest
from ase.build import bulk


class _StubRootstockCalculator:
    """Minimal stand-in for rootstock.RootstockCalculator."""

    entered = False
    exited = False
    call_kwargs: dict = {}

    def __init__(self, **kwargs):
        _StubRootstockCalculator.call_kwargs = kwargs
        _StubRootstockCalculator.entered = False
        _StubRootstockCalculator.exited = False
        self.parameters = {}
        n = kwargs.get("_n_atoms", 1)
        self.results = {
            "energy": -4.0,
            "forces": np.zeros((n, 3)),
            "stress": np.zeros(6),
        }
        self.implemented_properties = ["energy", "forces", "stress"]

    def __enter__(self):
        _StubRootstockCalculator.entered = True
        return self

    def __exit__(self, *args):
        _StubRootstockCalculator.exited = True
        return False

    def get_potential_energy(self):
        return self.results["energy"]

    def get_forces(self):
        return self.results["forces"]

    def get_stress(self):
        return self.results["stress"]

    def calculate(self, atoms, properties, system_changes):
        pass


def _make_rootstock_module() -> types.ModuleType:
    mod = types.ModuleType("rootstock")
    mod.__version__ = "0.0.0"
    mod.RootstockCalculator = _StubRootstockCalculator
    return mod


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _clear_cache():
    """Clear pick_calculator's lru_cache before and after every test."""
    from quacc.recipes.mlp._base import pick_calculator

    pick_calculator.__wrapped__.cache_clear()
    yield
    pick_calculator.__wrapped__.cache_clear()


@pytest.fixture
def rootstock_module(monkeypatch):
    """Inject a fake rootstock module so no real package is needed."""
    mod = _make_rootstock_module()
    monkeypatch.setitem(sys.modules, "rootstock", mod)
    return mod


def test_pick_calculator_rootstock(rootstock_module):
    """pick_calculator returns an entered RootstockCalculator."""
    from quacc.recipes.mlp._base import pick_calculator

    calc = pick_calculator(
        "rootstock", cluster="della", checkpoint="mace-mp-0-medium", device="cuda"
    )

    assert _StubRootstockCalculator.call_kwargs == {
        "cluster": "della",
        "checkpoint": "mace-mp-0-medium",
        "device": "cuda",
    }
    assert _StubRootstockCalculator.entered, "__enter__ must be called"
    assert isinstance(calc, _StubRootstockCalculator)


def test_pick_calculator_rootstock_cached(rootstock_module, monkeypatch):
    """The same kwargs return the cached object — constructor runs only once."""
    call_count = [0]
    original_init = _StubRootstockCalculator.__init__

    def counting_init(self, **kwargs):
        call_count[0] += 1
        original_init(self, **kwargs)

    monkeypatch.setattr(_StubRootstockCalculator, "__init__", counting_init)

    from quacc.recipes.mlp._base import pick_calculator

    pick_calculator("rootstock", cluster="della", checkpoint="mace-mp-0-medium")
    pick_calculator("rootstock", cluster="della", checkpoint="mace-mp-0-medium")

    assert call_count[0] == 1


def test_pick_calculator_rootstock_setup_kwargs(rootstock_module):
    """setup_kwargs are forwarded as a top-level kwarg to RootstockCalculator."""
    from quacc.recipes.mlp._base import pick_calculator

    pick_calculator(
        "rootstock",
        cluster="perlmutter",
        checkpoint="uma-s-1p1",
        setup_kwargs={"task": "omol"},
    )

    assert _StubRootstockCalculator.call_kwargs.get("setup_kwargs") == {"task": "omol"}


def test_static_job_rootstock(tmp_path, monkeypatch, rootstock_module):
    """static_job runs end-to-end with a stubbed RootstockCalculator."""
    monkeypatch.chdir(tmp_path)

    from quacc.recipes.mlp.core import static_job

    atoms = bulk("Cu")
    output = static_job(
        atoms,
        method="rootstock",
        cluster="della",
        checkpoint="mace-mp-0-medium",
        device="cuda",
    )

    assert "results" in output
    assert "atoms" in output
    assert output["atoms"] == atoms


def test_pick_calculator_rootstock_missing_package(monkeypatch):
    """pick_calculator raises ImportError when rootstock is not installed."""
    from quacc.recipes.mlp._base import pick_calculator

    monkeypatch.setitem(sys.modules, "rootstock", None)

    with pytest.raises((ImportError, ModuleNotFoundError)):
        pick_calculator("rootstock", cluster="della", checkpoint="mace-mp-0-medium")


def test_pick_calculator_unknown_method():
    """pick_calculator raises ValueError for an unrecognised method string."""
    from quacc.recipes.mlp._base import pick_calculator

    with pytest.raises(ValueError, match="Unrecognized"):
        pick_calculator("nonexistent-mlip")
