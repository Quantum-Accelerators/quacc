"""Tests for the MatterSim-compatible VASP recipes."""

from __future__ import annotations

from ase.build import bulk

from quacc.recipes.vasp.mattersim import mattersim_static_job


def test_mattersim_static_job(patch_nonmetallic_taskdoc):
    atoms = bulk("Ni") * (2, 1, 1)
    atoms[0].symbol = "O"
    del atoms.arrays["initial_magmoms"]

    output = mattersim_static_job(atoms, ncore=1)

    assert output["parameters"] == {
        "algo": "fast",
        "ediff": 0.0001,
        "encut": 520,
        "gamma": True,
        "isif": 3,
        "ismear": -5,
        "ispin": 2,
        "kpts": (5, 11, 11),
        "lasph": True,
        "ldau": True,
        "ldauj": [0, 0],
        "ldaul": [0, 2],
        "ldauprint": 1,
        "ldautype": 2,
        "ldauu": [0, 6.2],
        "lmaxmix": 4,
        "lorbit": 11,
        "lreal": "auto",
        "lwave": False,
        "magmom": [0.6, 5],
        "nelm": 100,
        "ncore": 1,
        "nsw": 0,
        "pp": "pbe",
        "pp_version": "original",
        "prec": "accurate",
        "setups": {"Ni": "_pv", "O": ""},
        "sigma": 0.05,
    }
    assert output["name"] == "MatterSim Static"
