"""Opt-in OC25 checkpoint integration tests; see the MLIP recipe guide."""

from __future__ import annotations

import gc
import os

import numpy as np
import pytest
from ase.build import add_adsorbate, fcc111, molecule
from ase.constraints import FixAtoms
from ase.io import read
from ase.optimize import FIRE

from quacc import change_settings
from quacc.recipes.mlip._base import pick_calculator
from quacc.recipes.mlip.core import relax_job, static_job
from quacc.recipes.mlip.md import md_job

pytestmark = pytest.mark.skipif(
    os.environ.get("QUACC_RUN_OC25_INTEGRATION") != "1",
    reason="Set QUACC_RUN_OC25_INTEGRATION=1 to run authenticated OC25 model tests.",
)


@pytest.fixture(scope="module", params=["esen", "uma"])
def model_kwargs(request):
    """Load each model only when explicitly requested, on the selected device."""
    torch = pytest.importorskip("torch")
    pytest.importorskip("fairchem.core")
    device = os.environ.get(
        "QUACC_OC25_DEVICE", "cuda" if torch.cuda.is_available() else "cpu"
    )
    if device not in {"cpu", "cuda"}:
        pytest.fail("QUACC_OC25_DEVICE must be 'cpu' or 'cuda'.")
    if device == "cuda" and not torch.cuda.is_available():
        pytest.fail("QUACC_OC25_DEVICE=cuda was requested, but CUDA is unavailable.")

    if request.param == "esen":
        kwargs = {
            "name_or_path": os.environ.get(
                "QUACC_OC25_ESEN_CHECKPOINT", "esen-sm-conserving-all-oc25"
            )
        }
    else:
        kwargs = {
            "name_or_path": os.environ.get("QUACC_OC25_UMA_CHECKPOINT", "uma-s-1p2"),
            "task_name": "oc25",
        }
    kwargs["device"] = device
    # Keep the smoke tests small and avoid compilation on either device. Use
    # the same inference settings for the independent and quacc calculations.
    kwargs["inference_settings"] = "batch"
    pick_calculator.__wrapped__.cache_clear()
    with change_settings({"FAIRCHEM_RAY_SERVE_BATCHING": False}):
        yield kwargs
    pick_calculator.__wrapped__.cache_clear()
    gc.collect()
    if device == "cuda":
        torch.cuda.empty_cache()


@pytest.fixture
def atoms():
    """Preserve an optional source structure or construct a small Cu/water slab."""
    if structure := os.environ.get("QUACC_OC25_STRUCTURE"):
        structure_atoms = read(structure)
    else:
        structure_atoms = fcc111("Cu", size=(3, 3, 3), vacuum=8.0)
        bottom = structure_atoms.get_tags() == 3
        structure_atoms.set_constraint(FixAtoms(mask=bottom))
        structure_atoms.set_tags(np.where(bottom, 0, 1))
        water = molecule("H2O")
        water.set_tags(2)
        add_adsorbate(structure_atoms, water, height=2.5, position="ontop")
        structure_atoms.pbc = True

    # A supplied structure without FixAtoms still needs a frozen atom to exercise
    # constraint handling; retain any other constraints the source provided.
    if not any(isinstance(c, FixAtoms) for c in structure_atoms.constraints):
        structure_atoms.set_constraint(
            [
                *structure_atoms.constraints,
                FixAtoms(indices=[int(np.argmin(structure_atoms.positions[:, 2]))]),
            ]
        )
    structure_atoms.calc = None
    return structure_atoms


def _fixed_indices(atoms):
    return np.unique(
        np.concatenate(
            [c.get_indices() for c in atoms.constraints if isinstance(c, FixAtoms)]
        )
    )


def _assert_preserved_structure(original, result):
    np.testing.assert_allclose(result.cell, original.cell, atol=1e-12, rtol=0)
    np.testing.assert_array_equal(result.pbc, original.pbc)
    np.testing.assert_array_equal(result.numbers, original.numbers)
    np.testing.assert_array_equal(result.get_tags(), original.get_tags())
    fixed = _fixed_indices(original)
    np.testing.assert_array_equal(_fixed_indices(result), fixed)
    np.testing.assert_allclose(
        result.positions[fixed], original.positions[fixed], atol=1e-12, rtol=0
    )


def test_static_matches_direct_fairchem(atoms, model_kwargs, tmp_path, monkeypatch):
    """Compare raw forces, including those on frozen atoms, and total energies."""
    from fairchem.core import FAIRChemCalculator

    monkeypatch.chdir(tmp_path)
    original = atoms.copy()
    direct = atoms.copy()
    direct.calc = FAIRChemCalculator.from_model_checkpoint(**model_kwargs)
    reference_energy = direct.get_potential_energy()
    reference_forces = direct.get_forces(apply_constraint=False).copy()
    # Release the independent reference calculator before quacc loads its model.
    direct.calc = None
    gc.collect()

    output = static_job(atoms, library="fairchem", **model_kwargs)
    assert np.isfinite(output["results"]["energy"])
    assert np.isfinite(output["results"]["forces"]).all()
    assert output["results"]["energy"] == pytest.approx(
        reference_energy, abs=1e-4, rel=0
    )
    np.testing.assert_allclose(
        output["results"]["forces"], reference_forces, atol=1e-4, rtol=0
    )
    _assert_preserved_structure(original, output["atoms"])
    np.testing.assert_array_equal(atoms.positions, original.positions)


def test_constrained_relaxation(atoms, model_kwargs, tmp_path, monkeypatch):
    """Exercise a bounded optimization; convergence is not asserted."""
    monkeypatch.chdir(tmp_path)
    original = atoms.copy()
    with change_settings({"CHECK_CONVERGENCE": False}):
        output = relax_job(
            atoms,
            library="fairchem",
            relax_cell=False,
            opt_params={"optimizer": FIRE, "max_steps": 3, "fmax": 0.05},
            **model_kwargs,
        )

    assert 1 <= len(output["trajectory"]) <= 4
    assert np.isfinite(output["results"]["energy"])
    assert np.isfinite(output["results"]["forces"]).all()
    for frame in output["trajectory"]:
        _assert_preserved_structure(original, frame)
    _assert_preserved_structure(original, output["atoms"])
    np.testing.assert_array_equal(atoms.positions, original.positions)


def test_fixed_cell_md(atoms, model_kwargs, tmp_path, monkeypatch):
    """Exercise a short, seeded trajectory without requesting stress."""
    monkeypatch.chdir(tmp_path)
    original = atoms.copy()
    output = md_job(
        atoms,
        library="fairchem",
        dynamics="nvt_langevin",
        steps=3,
        timestep_fs=0.1,
        temperature_K=100,
        md_params={
            "maxwell_boltzmann_kwargs": {"rng": np.random.default_rng(42)},
            "dynamics_kwargs": {"rng": np.random.default_rng(43)},
        },
        **model_kwargs,
    )

    assert len(output["trajectory"]) == 4
    assert np.isfinite(output["results"]["energy"])
    assert np.isfinite(output["results"]["forces"]).all()
    for frame in output["trajectory"]:
        _assert_preserved_structure(original, frame)
        assert np.isfinite(frame.get_potential_energy())
    _assert_preserved_structure(original, output["atoms"])
    np.testing.assert_array_equal(atoms.positions, original.positions)
    np.testing.assert_array_equal(atoms.get_momenta(), original.get_momenta())
