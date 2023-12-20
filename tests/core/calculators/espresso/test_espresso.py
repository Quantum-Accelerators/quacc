from __future__ import annotations

from ase import Atoms

from quacc.calculators.espresso.espresso import Espresso


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


def test_espresso_presets():
    calc_defaults = {
        "input_data": {
            "system": {"ecutwfc": 20, "ecutrho": 80, "occupations": "fixed"},
            "electrons": {"scf_must_converge": False},
        }
    }

    input_data = {"system": {"ecutwfc": 30}, "electrons": {"conv_thr": 1.0e-16}}

    preset = "esm_metal_slab_efficiency"

    atoms = Atoms(symbols="LiLaOZr")

    calc = Espresso(
        input_atoms=atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        input_data=input_data,
    )

    expected_parameters = {
        "input_data": {
            "control": {},
            "system": {
                "ecutwfc": 30,
                "ecutrho": 400.0,
                "occupations": "smearing",
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
            "ions": {},
            "cell": {},
            "rism": {},
        },
        "kspacing": 0.2,
        "pseudopotentials": {
            "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
            "Zr": "Zr_pbe_v1.uspp.F.UPF",
            "Li": "li_pbe_v1.4.uspp.F.UPF",
            "La": "La.paw.z_11.atompaw.wentzcovitch.v1.2.upf",
        },
    }
    assert calc.template.binary == "pw"
    assert calc.parameters == expected_parameters
