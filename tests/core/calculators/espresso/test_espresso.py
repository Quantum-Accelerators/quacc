from __future__ import annotations

from pathlib import Path

import pytest
from ase.atoms import Atoms

from quacc.calculators.espresso.espresso import Espresso, EspressoTemplate


def test_espresso_kwargs_handler():
    kwargs = {
        "kpts": (1, 1, 1),
        "input_data": {"system": {"ecutwfc": 30, "ecutrho": 80}},
    }

    preset = "sssp_1.3.0_pbe_efficiency"

    atoms = Atoms(symbols="LiLaOZr")

    calc = Espresso(input_atoms=atoms, preset=preset, **kwargs)

    expected_parameters = {
        "input_data": {"system": {"ecutwfc": 30, "ecutrho": 80}},
        "kpts": (1, 1, 1),
        "pseudopotentials": {
            "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
            "Zr": "Zr_pbe_v1.uspp.F.UPF",
            "Li": "li_pbe_v1.4.uspp.F.UPF",
            "La": "La.paw.z_11.atompaw.wentzcovitch.v1.2.upf",
        },
    }

    assert calc.template.binary == "pw"
    assert calc.parameters["kpts"] == expected_parameters["kpts"]
    assert (
        calc.parameters["input_data"]["system"]
        == expected_parameters["input_data"]["system"]
    )


def test_espresso_kwargs_handler_v2():
    preset = "sssp_1.3.0_pbe_efficiency"

    atoms = Atoms(symbols="LiLaOZr")

    calc = Espresso(input_atoms=atoms, preset=preset, kpts=None)

    expected_parameters = {
        "input_data": {"system": {"ecutwfc": 50.0, "ecutrho": 400.0}},
        "kpts": None,
        "pseudopotentials": {
            "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
            "Zr": "Zr_pbe_v1.uspp.F.UPF",
            "Li": "li_pbe_v1.4.uspp.F.UPF",
            "La": "La.paw.z_11.atompaw.wentzcovitch.v1.2.upf",
        },
    }

    assert calc.template.binary == "pw"
    assert calc.parameters["kpts"] == expected_parameters["kpts"]
    assert (
        calc.parameters["input_data"]["system"]
        == expected_parameters["input_data"]["system"]
    )
    assert (
        calc.parameters["pseudopotentials"] == expected_parameters["pseudopotentials"]
    )


def test_espresso_presets():
    input_data = {
        "system": {"ecutwfc": 200, "occupations": "fixed"},
        "electrons": {"scf_must_converge": False, "conv_thr": 1.0e-16},
    }

    preset = "esm_metal_slab_efficiency"

    atoms = Atoms(symbols="LiLaOZr")
    atoms.set_cell([5, 2, 10])

    calc = Espresso(
        input_atoms=atoms, preset=preset, kpts=[7, 17, 1], input_data=input_data
    )

    expected_parameters = {
        "input_data": {
            "system": {
                "ecutwfc": 200,
                "ecutrho": 400.0,
                "occupations": "fixed",
                "smearing": "marzari-vanderbilt",
                "degauss": 0.01,
                "assume_isolated": "esm",
                "esm_bc": "bc1",
            },
            "electrons": {
                "scf_must_converge": False,
                "conv_thr": 1.0e-16,
                "mixing_mode": "local-TF",
                "mixing_beta": 0.35,
            },
        },
        "kpts": [7, 17, 1],
        "pseudopotentials": {
            "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
            "Zr": "Zr_pbe_v1.uspp.F.UPF",
            "Li": "li_pbe_v1.4.uspp.F.UPF",
            "La": "La.paw.z_11.atompaw.wentzcovitch.v1.2.upf",
        },
    }
    assert calc.template.binary == "pw"
    assert (
        calc.parameters["input_data"]["system"]
        == expected_parameters["input_data"]["system"]
    )
    assert (
        calc.parameters["input_data"]["electrons"]
        == expected_parameters["input_data"]["electrons"]
    )
    assert calc.parameters["kpts"] == expected_parameters["kpts"]
    assert (
        calc.parameters["pseudopotentials"] == expected_parameters["pseudopotentials"]
    )


def test_espresso_presets_gamma():
    preset = "molecule_efficiency"

    calc = Espresso(preset=preset)

    assert calc.parameters["kpts"] is None


def test_espresso_bad_kpts():
    with pytest.raises(ValueError, match="Cannot specify both kpts and kspacing"):
        Espresso(kpts=(1, 1, 1), kspacing=0.5)


def test_espresso_kpts():
    calc = Espresso(kspacing=0.001, preset="metal_precision")

    assert "kpts" not in calc.parameters
    assert calc.parameters["kspacing"] == 0.001


def test_output_handler(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    parameters = {"input_data": {"control": {"outdir": "../test3/test2/test1"}}}

    test_path = "../test3/test2/test1"
    fake_template = EspressoTemplate()
    new_parameters = fake_template._output_handler(parameters, Path())

    assert str(new_parameters["input_data"]["control"]["outdir"]) != test_path

    test_path = "/test3/test2/test1"
    parameters["input_data"]["control"]["outdir"] = test_path
    fake_template = EspressoTemplate()

    new_parameters = fake_template._output_handler(parameters, Path())

    assert str(new_parameters["input_data"]["control"]["outdir"]) != test_path


def test_bad_calculator_params():
    atoms = Atoms(symbols="LiLaOZr")

    with pytest.raises(NotImplementedError, match="does not support the directory"):
        Espresso(input_atoms=atoms, kpts=(1, 1, 1), directory="bad")
