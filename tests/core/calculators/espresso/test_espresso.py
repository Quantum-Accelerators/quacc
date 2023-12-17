from __future__ import annotations

from pathlib import Path

import pytest
from ase import Atoms

from quacc.calculators.espresso.espresso import Espresso, EspressoTemplate


def test_espresso_kwargs_handler():
    calc_defaults = {
        "input_data": {"system": {"ecutwfc": 20, "ecutrho": 80, "conv_thr": 1e-8}}
    }

    input_data = {"system": {"ecutwfc": 30}}

    preset = "sssp_1.3.0_pbe_efficiency"

    atoms = Atoms(symbols="LiLaOZr")

    calc = Espresso(
        input_atoms=atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        input_data=input_data,
        kpts=(1, 1, 1),
    )

    expected_parameters = {
        "input_data": {
            "control": {},
            "system": {"ecutwfc": 30, "ecutrho": 400, "conv_thr": 1e-8},
            "electrons": {},
            "ions": {},
            "cell": {},
            "rism": {},
        },
        "kpts": (1, 1, 1),
        "pseudopotentials": {
            "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
            "Zr": "Zr_pbe_v1.uspp.F.UPF",
            "Li": "li_pbe_v1.4.uspp.F.UPF",
            "La": "La.paw.z_11.atompaw.wentzcovitch.v1.2.upf",
        },
    }

    assert calc.template.binary == "pw"
    assert calc.parameters == expected_parameters


def test_outdir_handler(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    parameters = {
        "input_data": {
            "control": {},
            "system": {"outdir": "../test3/test2/test1"},
            "electrons": {},
            "ions": {},
            "cell": {},
            "rism": {},
        }
    }

    test_path = "../test3/test2/test1"
    fake_template = EspressoTemplate()
    new_parameters = fake_template._outdir_handler(parameters)
    assert len(fake_template.outdirs) == 1
    assert new_parameters["input_data"]["system"]["outdir"] == test_path
    assert not Path(test_path).exists()

    test_path = "/test3/test2/test1"
    parameters["input_data"]["system"]["outdir"] = test_path
    fake_template = EspressoTemplate()
    new_parameters = fake_template._outdir_handler(parameters)
    assert len(fake_template.outdirs) == 1
    assert new_parameters["input_data"]["system"]["outdir"] == test_path
    assert not Path(test_path).exists()

    test_path = Path("test3/test2/test1")
    parameters["input_data"]["system"]["outdir"] = test_path
    fake_template = EspressoTemplate()
    new_parameters = fake_template._outdir_handler(parameters)
    assert len(fake_template.outdirs) == 0
    assert (
        new_parameters["input_data"]["system"]["outdir"]
        == test_path.expanduser().resolve()
    )
    assert Path(tmp_path, "test3/test2/test1").exists()


def test_bad_calculator_params():
    calc_defaults = {
        "input_data": {"system": {"ecutwfc": 20, "ecutrho": 80, "conv_thr": 1e-8}}
    }

    input_data = {"system": {"ecutwfc": 30}}

    atoms = Atoms(symbols="LiLaOZr")

    with pytest.raises(ValueError):
        Espresso(
            input_atoms=atoms,
            calc_defaults=calc_defaults,
            input_data=input_data,
            kpts=(1, 1, 1),
            directory="bad",
        )
