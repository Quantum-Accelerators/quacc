"""Tests for the MatterSim-compatible VASP recipes."""

from __future__ import annotations

from ase.build import bulk

from quacc.recipes.vasp.mattersim import mattersim_static_job


def test_mattersim_static_job(monkeypatch):
    atoms = bulk("Al")

    monkeypatch.setattr(
        "quacc.recipes.vasp.mattersim.run_and_summarize",
        lambda atoms, **kwargs: {"atoms": atoms, **kwargs},
    )

    output = mattersim_static_job(atoms, encut=600)

    assert output["calc_defaults"] == {
        "algo": "Fast",
        "ediff": 5e-05,
        "encut": 520.0,
        "gamma": True,
        "ibrion": 2,
        "isif": 3,
        "ismear": -5,
        "ispin": 2,
        "kpts": (9, 9, 9),
        "lasph": True,
        "lorbit": 11,
        "lreal": "Auto",
        "lwave": False,
        "magmom": [0.6],
        "nelm": 100,
        "nsw": 0,
        "pp": "PBE",
        "pp_version": "original",
        "prec": "Accurate",
        "setups": {"Al": ""},
        "sigma": 0.05,
    }
    assert output["calc_swaps"] == {"encut": 600}
    assert output["additional_fields"] == {"name": "MatterSim Static"}
