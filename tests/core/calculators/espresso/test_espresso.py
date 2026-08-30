from __future__ import annotations

from pathlib import Path

import pytest
from ase.atoms import Atoms

from quacc.calculators.espresso.espresso import Espresso, EspressoTemplate
from quacc.calculators.espresso.utils import grid_copy_files, prepare_copy_files


def test_grid_copy_files_non_gamma():
    directory = Path("previous")

    files = grid_copy_files({"inputph": {}}, directory, 1, (0.5, 0.0, 0.0))[directory]

    assert Path("_ph0", "pwscf.wfc*") in files


@pytest.mark.parametrize("binary", ["dos", "fs"])
def test_prepare_copy_files_density_only(binary):
    to_copy = prepare_copy_files({}, binary=binary)

    assert Path("pwscf.save", "charge-density.*") in to_copy


@pytest.mark.parametrize("binary", ["projwfc", "bands"])
def test_prepare_copy_files_with_wavefunctions(binary):
    to_copy = prepare_copy_files({}, binary=binary)

    assert Path("pwscf.save", "wfc*.*") in to_copy


def test_prepare_copy_files_pp_wavefunctions():
    to_copy = prepare_copy_files({"input_data": {"plot_num": 3}}, binary="pp")

    assert Path("pwscf.save", "wfc*.*") in to_copy


def test_prepare_copy_files_ph_interpolation():
    parameters = {
        "input_data": {"inputph": {"lqdir": True, "ldvscf_interpolate": True}}
    }

    to_copy = prepare_copy_files(parameters, binary="ph")

    assert Path("_ph*", "pwscf.q_*", "pwscf.dvscf*") in to_copy


def test_prepare_copy_files_postahc():
    to_copy = prepare_copy_files({}, binary="postahc")

    assert Path("pwscf.save", "data-file-schema.*") in to_copy


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


def test_search_keyword_missing():
    assert EspressoTemplate._search_keyword({"input_data": {}}, "prefix") is None


def test_ph_grid_recovery_directory(tmp_path):
    template = EspressoTemplate(binary="ph")
    parameters = {
        "input_data": {
            "inputph": {
                "qplot": True,
                "lqdir": True,
                "recover": True,
                "start_q": 2,
                "prefix": "test",
                "outdir": str(tmp_path),
            }
        }
    }

    template._sanity_checks(parameters)

    assert (tmp_path / "_ph0" / "test.q_1").is_dir()


def test_espresso_preset_without_pseudopotentials(tmp_path):
    preset = tmp_path / "preset.yaml"
    preset.write_text("input_data:\n  system:\n    occupations: fixed\n")

    calc = Espresso(preset=preset, kpts=None)

    assert calc.parameters["input_data"]["system"]["occupations"] == "fixed"


def test_bad_calculator_params():
    atoms = Atoms(symbols="LiLaOZr")

    with pytest.raises(NotImplementedError, match="does not support the directory"):
        Espresso(input_atoms=atoms, kpts=(1, 1, 1), directory="bad")


def test_prepare_copy_files_pw_restart():
    parameters = {
        "input_data": {"control": {"restart_mode": "restart"}, "electrons": {}}
    }

    to_copy = prepare_copy_files(parameters)

    assert Path("pwscf.wfc*") in to_copy
    assert Path("pwscf.mix*") in to_copy
    assert Path("pwscf.restart_k*") in to_copy
    assert Path("pwscf.restart_scf*") in to_copy
    assert Path("pwscf.save", "charge-density.*") in to_copy
    assert Path("pwscf.save", "wfc*.*") in to_copy


def test_prepare_copy_files_ph_recovery():
    parameters = {"input_data": {"inputph": {"lqdir": True, "recover": True}}}

    to_copy = prepare_copy_files(parameters, binary="ph")

    assert Path("_ph*", "pwscf.q_*", "pwscf.save", "charge-density.*") in to_copy
    assert Path("_ph*", "pwscf.phsave") in to_copy


@pytest.mark.parametrize(
    ("binary", "expected"),
    [
        ("matdyn", Path("q2r.fc*")),
        ("q2r", Path("matdyn*")),
        ("dvscf_q2r", Path("matdyn0*")),
    ],
)
def test_prepare_copy_files_postprocessing(binary, expected):
    assert expected in prepare_copy_files({}, binary=binary)
