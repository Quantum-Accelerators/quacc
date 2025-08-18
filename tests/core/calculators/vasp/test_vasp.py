from __future__ import annotations

import os
from copy import deepcopy
from logging import INFO, WARNING, getLogger
from pathlib import Path
from shutil import which

import numpy as np
import pytest
from ase.atoms import Atoms
from ase.build import bulk
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.calculators.vasp import Vasp as Vasp_
from ase.constraints import FixAtoms, FixBondLength
from ase.io import read
from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet

from quacc import change_settings, get_settings
from quacc.calculators.vasp import Vasp, presets
from quacc.calculators.vasp.params import MPtoASEConverter
from quacc.schemas.prep import prep_next_run

FILE_DIR = Path(__file__).parent
PSEUDO_DIR = FILE_DIR / "fake_pseudos"
LOGGER = getLogger(__name__)
LOGGER.propagate = True


@pytest.fixture
def atoms_mag():
    return read(FILE_DIR / "OUTCAR_mag.gz")


@pytest.fixture
def atoms_nomag():
    return read(FILE_DIR / "OUTCAR_nomag.gz")


@pytest.fixture
def atoms_nospin():
    return read(FILE_DIR / "OUTCAR_nospin.gz")


def test_vanilla_vasp():
    atoms = bulk("Cu")
    calc = Vasp(atoms, incar_copilot=False)
    assert calc.asdict() == Vasp_().asdict()

    atoms = bulk("Cu")
    calc = Vasp(atoms, use_custodian=False, incar_copilot=False)
    assert calc.asdict() == Vasp_().asdict()

    atoms = bulk("Cu")
    calc = Vasp(atoms, use_custodian=False, kspacing=0.5, incar_copilot=False)
    assert calc.input_params["gamma"] is None
    assert calc.kpts is None

    atoms = bulk("Cu")
    calc = Vasp(atoms, encut=None, incar_copilot=False)
    assert calc.asdict() == Vasp_().asdict()


@pytest.mark.parametrize(
    "preset",
    [
        "BulkSet",
        Path(presets.__file__).parent / "BulkSet.yaml",
        str(Path(presets.__file__).parent / "BulkSet.yaml"),
    ],
)
def test_presets_basic(preset):
    atoms = bulk("Co") * (2, 2, 1)
    atoms[-1].symbol = "Fe"

    calc = Vasp(atoms, preset=preset)
    atoms.calc = calc
    assert calc.xc.lower() == "pbe"
    assert calc.string_params["algo"] == "fast"
    assert calc.exp_params["ediff"] == 1e-5
    assert calc.float_params["encut"] == 520


def test_presets2():
    atoms = bulk("Co") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    calc = Vasp(atoms, xc="rpbe", preset="SlabSet")
    assert calc.xc.lower() == "rpbe"
    assert calc.string_params["algo"] == "fast"
    assert calc.exp_params["ediff"] == 1e-5
    assert calc.float_params["encut"] == 450


def test_presets_mp():
    atoms = bulk("Co") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    parameters = MPtoASEConverter(atoms=atoms).convert_dict_set(MPScanRelaxSet)
    calc = Vasp(atoms, xc="scan", **parameters)
    assert calc.xc.lower() == "scan"
    assert calc.string_params["algo"].lower() == "all"
    assert calc.exp_params["ediff"] == 1e-5


def test_rosen_preset1():
    params = {
        "algo": "all",
        "ediff": 1e-06,
        "ediffg": -0.02,
        "efermi": "midgap",
        "encut": 520,
        "gamma": None,
        "gga": "PE",
        "gga_compat": False,
        "ibrion": 2,
        "ismear": 0,
        "ivdw": 12,
        "kpts": None,
        "kspacing": 0.2,
        "lasph": True,
        "lcharg": False,
        "lmaxmix": 4,
        "lorbit": 11,
        "lreal": False,
        "lwave": False,
        "nelm": 150,
        "nelmin": 3,
        "nsw": 100,
        "pp": "PBE",
        "prec": "accurate",
        "setups": {
            "Ac": "",
            "Ag": "",
            "Al": "",
            "Am": "",
            "Ar": "",
            "As": "",
            "At": "",
            "Au": "",
            "B": "",
            "Ba": "_sv",
            "Be": "",
            "Bi": "_d",
            "Br": "",
            "C": "",
            "Ca": "_sv",
            "Cd": "",
            "Ce": "",
            "Cf": "",
            "Cl": "",
            "Cm": "",
            "Co": "",
            "Cr": "_pv",
            "Cs": "_sv",
            "Cu": "",
            "Dy": "_3",
            "Er": "_3",
            "Eu": "_3",
            "F": "",
            "Fe": "",
            "Fr": "_sv",
            "Ga": "_d",
            "Gd": "_3",
            "Ge": "_d",
            "H": "",
            "He": "",
            "Hf": "_pv",
            "Hg": "",
            "Ho": "_3",
            "I": "",
            "In": "_d",
            "Ir": "",
            "K": "_sv",
            "Kr": "",
            "La": "",
            "Li": "_sv",
            "Lu": "_3",
            "Mg": "",
            "Mn": "_pv",
            "Mo": "_sv",
            "N": "",
            "Na": "_pv",
            "Nb": "_sv",
            "Nd": "_3",
            "Ne": "",
            "Ni": "",
            "Np": "",
            "O": "",
            "Os": "",
            "P": "",
            "Pa": "",
            "Pb": "_d",
            "Pd": "",
            "Pm": "_3",
            "Po": "_d",
            "Pr": "_3",
            "Pt": "",
            "Pu": "",
            "Ra": "_sv",
            "Rb": "_sv",
            "Re": "",
            "Rh": "_pv",
            "Rn": "",
            "Ru": "_pv",
            "S": "",
            "Sb": "",
            "Sc": "_sv",
            "Se": "",
            "Si": "",
            "Sm": "_3",
            "Sn": "_d",
            "Sr": "_sv",
            "Ta": "_pv",
            "Tb": "_3",
            "Tc": "_pv",
            "Te": "",
            "Th": "",
            "Ti": "_sv",
            "Tl": "_d",
            "Tm": "_3",
            "U": "",
            "V": "_sv",
            "W": "_sv",
            "Xe": "",
            "Y": "_sv",
            "Yb": "_3",
            "Zn": "",
            "Zr": "_sv",
        },
        "sigma": 0.05,
        "xc": "pbe",
    }

    atoms = bulk("Cu")

    calc = Vasp(atoms, preset="RosenSetPBE", nsw=100)
    calc.parameters.pop("ncore")
    assert calc.parameters == params

    calc = Vasp(atoms, preset="RosenFastSetPBE", nsw=100)
    calc.parameters.pop("ncore")
    params2 = deepcopy(params)
    params2.pop("kspacing")
    params2["kpts"] = [12, 12, 12]
    params2["gamma"] = True
    params2["ediff"] = 1e-5
    assert calc.parameters == params2


def test_rosen_preset2():
    params = {
        "algo": "all",
        "ediff": 1e-06,
        "ediffg": -0.02,
        "efermi": "midgap",
        "encut": 520,
        "gamma": None,
        "gga_compat": False,
        "ibrion": 2,
        "ismear": 0,
        "ivdw": 13,
        "kpts": None,
        "kspacing": 0.2,
        "lasph": True,
        "lcharg": False,
        "lmaxmix": 4,
        "lorbit": 11,
        "lreal": False,
        "lwave": False,
        "metagga": "R2SCAN",
        "nelm": 150,
        "nelmin": 3,
        "nsw": 100,
        "pp": "PBE",
        "prec": "accurate",
        "setups": {
            "Ac": "",
            "Ag": "",
            "Al": "",
            "Am": "",
            "Ar": "",
            "As": "",
            "At": "",
            "Au": "",
            "B": "",
            "Ba": "_sv",
            "Be": "",
            "Bi": "_d",
            "Br": "",
            "C": "",
            "Ca": "_sv",
            "Cd": "",
            "Ce": "",
            "Cf": "",
            "Cl": "",
            "Cm": "",
            "Co": "",
            "Cr": "_pv",
            "Cs": "_sv",
            "Cu": "",
            "Dy": "_3",
            "Er": "_3",
            "Eu": "_3",
            "F": "",
            "Fe": "",
            "Fr": "_sv",
            "Ga": "_d",
            "Gd": "_3",
            "Ge": "_d",
            "H": "",
            "He": "",
            "Hf": "_pv",
            "Hg": "",
            "Ho": "_3",
            "I": "",
            "In": "_d",
            "Ir": "",
            "K": "_sv",
            "Kr": "",
            "La": "",
            "Li": "_sv",
            "Lu": "_3",
            "Mg": "",
            "Mn": "_pv",
            "Mo": "_sv",
            "N": "",
            "Na": "_pv",
            "Nb": "_sv",
            "Nd": "_3",
            "Ne": "",
            "Ni": "",
            "Np": "",
            "O": "",
            "Os": "",
            "P": "",
            "Pa": "",
            "Pb": "_d",
            "Pd": "",
            "Pm": "_3",
            "Po": "_d",
            "Pr": "_3",
            "Pt": "",
            "Pu": "",
            "Ra": "_sv",
            "Rb": "_sv",
            "Re": "",
            "Rh": "_pv",
            "Rn": "",
            "Ru": "_pv",
            "S": "",
            "Sb": "",
            "Sc": "_sv",
            "Se": "",
            "Si": "",
            "Sm": "_3",
            "Sn": "_d",
            "Sr": "_sv",
            "Ta": "_pv",
            "Tb": "_3",
            "Tc": "_pv",
            "Te": "",
            "Th": "",
            "Ti": "_sv",
            "Tl": "_d",
            "Tm": "_3",
            "U": "",
            "V": "_sv",
            "W": "_sv",
            "Xe": "",
            "Y": "_sv",
            "Yb": "_3",
            "Zn": "",
            "Zr": "_sv",
        },
        "sigma": 0.05,
        "vdw_a1": 0.51559235,
        "vdw_a2": 5.77342911,
        "vdw_s6": 1.0,
        "vdw_s8": 0.6018749,
        "xc": "r2scan",
    }
    atoms = bulk("Cu")
    calc = Vasp(atoms, preset="RosenSetR2SCAN", nsw=100)
    calc.parameters.pop("ncore")
    assert calc.parameters == params

    calc = Vasp(atoms, preset="RosenFastSetR2SCAN", nsw=100)
    calc.parameters.pop("ncore")
    params2 = deepcopy(params)
    params2.pop("kspacing")
    params2["kpts"] = [12, 12, 12]
    params2["gamma"] = True
    params2["ediff"] = 1e-5
    assert calc.parameters == params2


def test_lmaxmix():
    atoms = bulk("Cu")
    calc = Vasp(atoms)
    assert calc.int_params["lmaxmix"] == 4

    atoms = bulk("Ce")
    calc = Vasp(atoms)
    assert calc.int_params["lmaxmix"] == 6

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[-1].symbol = "Ce"
    calc = Vasp(atoms)
    assert calc.int_params["lmaxmix"] == 6


def test_li_sv(caplog):
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].symbol = "Li"
    with caplog.at_level(WARNING):
        Vasp(atoms, isif=3, nsw=2, encut=1000, setups={"Li": "Li_sv"})
    assert "Pulay" not in caplog.text


def test_li_sv2(caplog):
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].symbol = "Li"
    with caplog.at_level(WARNING):
        Vasp(atoms, isif=3, nsw=2, encut=400, setups={"Li": "Li_sv"})
    assert "Pulay" in caplog.text


def test_li_sv3(caplog):
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].symbol = "Li"
    with caplog.at_level(WARNING):
        Vasp(atoms, isif=3, nsw=2, setups={"Li": "Li_sv"})
    assert "Pulay" in caplog.text


def test_li_sv4(caplog):
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].symbol = "Li"
    with caplog.at_level(WARNING):
        Vasp(atoms, isif=3, nsw=2, preset="BulkSet")
    assert "Pulay" in caplog.text


def test_li_sv5(caplog):
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].symbol = "Li"
    with caplog.at_level(WARNING):
        Vasp(atoms, isif=3, nsw=2, setups={"Li": ""})
    assert "Pulay" in caplog.text


def test_he(caplog):
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].symbol = "He"
    with caplog.at_level(WARNING):
        Vasp(atoms, isif=3, nsw=2, encut=1000)
    assert "Pulay" not in caplog.text


def test_he2(caplog):
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].symbol = "He"
    with caplog.at_level(WARNING):
        Vasp(atoms, isif=3, nsw=2, encut=400)
    assert "Pulay" in caplog.text


def test_he3(caplog):
    atoms = bulk("Cu") * (2, 2, 2)
    atoms[0].symbol = "He"
    with caplog.at_level(WARNING):
        Vasp(atoms, isif=3, nsw=2)
    assert "Pulay" in caplog.text


def test_no_encut(caplog):
    atoms = bulk("Cu")
    with caplog.at_level(WARNING):
        Vasp(atoms)
    assert "Pulay" not in caplog.text


def test_no_encut2(caplog):
    atoms = bulk("Cu")

    with caplog.at_level(WARNING):
        Vasp(atoms, isif=3, nsw=2)
    assert "Pulay" in caplog.text


def test_no_encut3(caplog):
    atoms = bulk("Cu")

    with caplog.at_level(WARNING):
        Vasp(atoms, isif=3, nsw=2, encut=30000)
    assert "Pulay" not in caplog.text


def test_autodipole():
    atoms = bulk("Cu")
    com = atoms.get_center_of_mass(scaled=True)
    calc = Vasp(atoms, auto_dipole=True)
    assert calc.bool_params["ldipol"] is True
    assert calc.int_params["idipol"] == 3
    assert np.array_equal(calc.list_float_params["dipol"], com)

    calc = Vasp(atoms, auto_dipole=True, idipol=2)
    assert calc.bool_params["ldipol"] is True
    assert calc.int_params["idipol"] == 2
    assert np.array_equal(calc.list_float_params["dipol"], com)

    calc = Vasp(atoms, preset="SlabSet")
    assert calc.bool_params["ldipol"] is True
    assert calc.int_params["idipol"] == 3
    assert np.array_equal(calc.list_float_params["dipol"], com)

    calc = Vasp(atoms, preset="SlabSet", idipol=2)
    assert calc.bool_params["ldipol"] is True
    assert calc.int_params["idipol"] == 2
    assert np.array_equal(calc.list_float_params["dipol"], com)

    calc = Vasp(atoms, auto_dipole=False, preset="SlabSet", idipol=2)
    assert calc.bool_params["ldipol"] is None
    assert calc.int_params["idipol"] == 2
    assert calc.list_float_params["dipol"] is None

    calc = Vasp(atoms, auto_dipole=False, preset="SlabSet", idipol=2)
    assert calc.bool_params["ldipol"] is None
    assert calc.int_params["idipol"] == 2
    assert calc.list_float_params["dipol"] is None

    calc = Vasp(atoms, lsorbit=True)
    assert calc.bool_params["lsorbit"] is True
    assert calc.int_params["isym"] == -1


def test_kspacing():
    atoms = bulk("Cu")
    calc = Vasp(atoms, kspacing=0.1, ismear=-5)
    assert calc.int_params["ismear"] == -5

    calc = Vasp(atoms, kspacing=100, ismear=-5)
    assert calc.int_params["ismear"] == -5

    calc = Vasp(atoms, kspacing=0.1, preset="BulkSet")
    assert calc.float_params["kspacing"] == 0.1
    assert calc.kpts is None

    calc = Vasp(atoms, kspacing=0.1, gamma=True)
    assert calc.float_params["kspacing"] == 0.1
    assert calc.input_params["gamma"] is None
    assert calc.kpts is None

    calc = Vasp(atoms, kpts=[1, 1, 1], kgamma=True)
    assert calc.bool_params["kgamma"] is None


def test_kspacing_aggressive():
    with change_settings({"VASP_INCAR_COPILOT": "aggressive"}):
        atoms = bulk("Cu")
        calc = Vasp(atoms, kspacing=0.1, ismear=-5)
        assert calc.int_params["ismear"] == -5

        calc = Vasp(atoms, kspacing=100, ismear=-5)
        assert calc.int_params["ismear"] == 0


def test_magmoms(atoms_mag, atoms_nomag, atoms_nospin):
    atoms = bulk("Mg")
    calc = Vasp(atoms)
    atoms.calc = calc
    assert atoms.has("initial_magmoms") is False

    atoms = bulk("Mg")
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    calc = Vasp(atoms)
    atoms.calc = calc
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * len(atoms)

    atoms = bulk("Cu") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert atoms.get_chemical_symbols() == ["Cu", "Cu", "Cu", "Fe"]
    assert atoms.get_initial_magnetic_moments().tolist() == [2.0] * (len(atoms) - 1) + [
        5.0
    ]

    atoms = bulk("Zn") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    calc = Vasp(atoms, preset="BulkSet", preset_mag_default=2.5)
    atoms.calc = calc
    assert atoms.get_initial_magnetic_moments().tolist() == [2.5] * (len(atoms) - 1) + [
        5.0
    ]

    atoms = bulk("Eu") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    calc = Vasp(atoms, preset="SlabSet")
    atoms.calc = calc
    assert atoms.get_initial_magnetic_moments().tolist() == [7.0] * (len(atoms) - 1) + [
        5.0
    ]

    atoms = bulk("Cu") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    parameters = MPtoASEConverter(atoms=atoms).convert_dict_set(MPScanRelaxSet)
    calc = Vasp(atoms, **parameters)
    atoms.calc = calc
    assert atoms.get_chemical_symbols() == ["Cu", "Cu", "Cu", "Fe"]
    assert calc.parameters["magmom"] == [0.6, 0.6, 0.6, 5.0]

    atoms = bulk("Cu") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    atoms.set_initial_magnetic_moments([3.14] * (len(atoms) - 1) + [1.0])
    parameters = MPtoASEConverter(atoms=atoms).convert_dict_set(MPScanRelaxSet)
    calc = Vasp(atoms, **parameters)
    atoms.calc = calc
    assert atoms.get_initial_magnetic_moments().tolist() == [3.14] * (
        len(atoms) - 1
    ) + [1.0]

    atoms = bulk("Co") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    atoms.set_initial_magnetic_moments([0.0] * len(atoms))
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 0)

    atoms = deepcopy(atoms_mag)
    mags = atoms.get_magnetic_moments()
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.array_equal(atoms.get_initial_magnetic_moments(), mags) is True

    atoms = deepcopy(atoms_mag)
    calc = Vasp(atoms, preset="BulkSet", mag_cutoff=2.0)
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 0)

    atoms = deepcopy(atoms_mag)
    assert atoms.get_magnetic_moments()[0] == 0.468
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert atoms.get_initial_magnetic_moments()[0] == 0.468

    atoms = deepcopy(atoms_nomag)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 0)
    atoms.calc.results = {"energy": -1.0}  # mock calculation run
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert atoms.has("initial_magmoms") is True
    assert np.all(atoms.get_initial_magnetic_moments() == 0)

    atoms = deepcopy(atoms_nomag)
    calc = SinglePointDFTCalculator(atoms, magmoms=[0.0] * len(atoms))
    atoms.calc = calc
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 0)

    atoms = deepcopy(atoms_nomag)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert atoms.has("initial_magmoms") is True
    assert np.all(atoms.get_initial_magnetic_moments() == 0)
    atoms.calc.results = {"magmoms": [0.0] * len(atoms)}
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 0)

    atoms = deepcopy(atoms_nomag)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert atoms.has("initial_magmoms") is True
    assert np.all(atoms.get_initial_magnetic_moments() == 0)
    atoms.calc.results = {"magmoms": [1.0] * len(atoms)}
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 1)

    atoms = deepcopy(atoms_nospin)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert atoms.has("initial_magmoms") is True
    assert np.all(atoms.get_initial_magnetic_moments() == 0)
    atoms.calc.results = {"magmoms": [1.0] * len(atoms)}
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 1)

    atoms = bulk("Cu")
    calc = Vasp(atoms, preset="BulkSet", copy_magmoms=False)
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 2.0)
    atoms.calc.results = {"magmoms": [3.0] * len(atoms)}
    calc = Vasp(atoms, preset="BulkSet", copy_magmoms=False)
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 2.0)

    atoms = bulk("Mg")
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 1.0)
    atoms.calc.results = {"magmoms": [0.0] * len(atoms)}
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 0)

    atoms = bulk("Mg")
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 1.0)
    atoms.calc.results = {"magmoms": [-0.01] * len(atoms)}
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 0)

    atoms = bulk("Mg")
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == 1.0)
    atoms.calc.results = {"magmoms": [-5] * len(atoms)}
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert np.all(atoms.get_initial_magnetic_moments() == -5)


def test_prep_magmoms2(atoms_mag, atoms_nomag, atoms_nospin):
    atoms = deepcopy(atoms_mag)
    mags = atoms.get_magnetic_moments()
    atoms = prep_next_run(atoms, move_magmoms=True)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert atoms.has("initial_magmoms") is True
    assert atoms.get_initial_magnetic_moments().tolist() == mags.tolist()

    atoms = deepcopy(atoms_nomag)
    atoms = prep_next_run(atoms, move_magmoms=True)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert atoms.has("initial_magmoms") is True
    assert np.all(atoms.get_initial_magnetic_moments() == 0)

    atoms = deepcopy(atoms_nospin)
    atoms = prep_next_run(atoms, move_magmoms=True)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    assert atoms.has("initial_magmoms") is True
    assert np.all(atoms.get_initial_magnetic_moments() == 0)

    atoms = deepcopy(atoms_mag)
    atoms = prep_next_run(atoms, move_magmoms=True)
    calc = Vasp(atoms, preset="BulkSet", mag_cutoff=10.0)
    atoms.calc = calc
    assert atoms.has("initial_magmoms") is True
    assert np.all(atoms.get_initial_magnetic_moments() == 0)

    atoms = bulk("Mg")
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"energy": -1.0, "magmoms": [0.0] * len(atoms)}
    atoms = prep_next_run(atoms, move_magmoms=True)
    atoms *= (2, 2, 2)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    np.all(atoms.get_initial_magnetic_moments() == 0.0)

    atoms = bulk("Mg")
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"energy": -1.0, "magmoms": [-0.02] * len(atoms)}
    atoms = prep_next_run(atoms, move_magmoms=True)
    atoms *= (2, 2, 2)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    np.all(atoms.get_initial_magnetic_moments() == 0.0)

    atoms = bulk("Mg")
    calc = Vasp(atoms)
    atoms.calc = calc
    atoms.calc.results = {"energy": -1.0}
    atoms = prep_next_run(atoms, move_magmoms=True)
    atoms *= (2, 2, 2)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    np.all(atoms.get_initial_magnetic_moments() == 0.0)

    atoms = bulk("Mg")
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    atoms.calc.results = {"energy": -1.0, "magmoms": [3.14] * len(atoms)}
    atoms = prep_next_run(atoms, move_magmoms=True)
    atoms *= (2, 2, 2)
    calc = Vasp(atoms, preset="BulkSet")
    atoms.calc = calc
    np.all(atoms.get_initial_magnetic_moments() == 3.14)


def test_unused_flags():
    atoms = bulk("Cu")
    calc = Vasp(atoms, preset="BulkSet", potim=1.5, nsw=0)
    assert calc.int_params["nsw"] == 0
    assert calc.exp_params["ediffg"] is None
    assert calc.int_params["isif"] is None
    assert calc.float_params["potim"] is None

    calc = Vasp(atoms, ldau=False, ldauprint=2)
    assert calc.int_params["ldauprint"] is None
    assert calc.bool_params["ldau"] is None


def test_lasph():
    atoms = bulk("Cu")

    calc = Vasp(atoms, xc="rpbe")
    assert calc.bool_params["lasph"] is None

    calc = Vasp(atoms, xc="m06l")
    assert calc.bool_params["lasph"] is True

    calc = Vasp(atoms, xc="m06l", lasph=False)
    assert calc.bool_params["lasph"] is False

    calc = Vasp(atoms, xc="hse06")
    assert calc.bool_params["lasph"] is True

    calc = Vasp(atoms, ldau_luj={"Cu": {"L": 2, "U": 5, "J": 0.0}})
    assert calc.bool_params["lasph"] is True


def test_lasph_aggressive():
    with change_settings({"VASP_INCAR_COPILOT": "aggressive"}):
        atoms = bulk("Cu")

        calc = Vasp(atoms, xc="rpbe")
        assert calc.bool_params["lasph"] is None

        calc = Vasp(atoms, xc="m06l")
        assert calc.bool_params["lasph"] is True

        calc = Vasp(atoms, xc="m06l", lasph=False)
        assert calc.bool_params["lasph"] is True

        calc = Vasp(atoms, xc="hse06")
        assert calc.bool_params["lasph"] is True

        calc = Vasp(atoms, ldau_luj={"Cu": {"L": 2, "U": 5, "J": 0.0}})
        assert calc.bool_params["lasph"] is True


def test_algo():
    atoms = bulk("Cu")

    calc = Vasp(atoms, xc="rpbe")
    assert calc.string_params["algo"] is None

    calc = Vasp(atoms, xc="m06l")
    assert calc.string_params["algo"] == "all"

    calc = Vasp(atoms, xc="hse06")
    assert calc.string_params["algo"] == "normal"


def test_algo_aggressive():
    with change_settings({"VASP_INCAR_COPILOT": "aggressive"}):
        atoms = bulk("Cu")

        calc = Vasp(atoms, xc="rpbe")
        assert calc.string_params["algo"] is None

        calc = Vasp(atoms, xc="m06l")
        assert calc.string_params["algo"] == "all"

        calc = Vasp(atoms, xc="m06l", algo="fast")
        assert calc.string_params["algo"] == "all"

        calc = Vasp(atoms, xc="hse06")
        assert calc.string_params["algo"] == "normal"


def test_kpar():
    atoms = bulk("Cu")

    calc = Vasp(atoms, kpts=[2, 2, 1], kpar=4)
    assert calc.int_params["kpar"] == 4

    with change_settings({"VASP_INCAR_COPILOT": "aggressive"}):
        calc = Vasp(atoms, kpts=[2, 2, 1], kpar=4)
        assert calc.int_params["kpar"] == 4

        calc = Vasp(atoms, kpar=4)
        assert calc.int_params["kpar"] == 1


def test_isym_aggressive():
    with change_settings({"VASP_INCAR_COPILOT": "aggressive"}):
        atoms = bulk("Cu")

        calc = Vasp(atoms, isym=2)
        assert calc.int_params["isym"] == 2

        calc = Vasp(atoms, isym=0)
        assert calc.int_params["isym"] == 0

        calc = Vasp(atoms, xc="hse06", isym=2)
        assert calc.int_params["isym"] == 3

        calc = Vasp(atoms, isym=2, nsw=100)
        assert calc.int_params["isym"] == 2

        calc = Vasp(atoms, xc="hse06", isym=2, nsw=100)
        assert calc.int_params["isym"] == 3


def test_ncore_aggressive():
    with change_settings({"VASP_INCAR_COPILOT": "aggressive"}):
        atoms = bulk("Cu")

        atoms *= (2, 2, 2)
        calc = Vasp(atoms, ncore=4)
        assert calc.int_params["ncore"] == 4

        calc = Vasp(atoms, ncore=4, lhfcalc=True)
        assert calc.int_params["ncore"] == 1

        calc = Vasp(atoms, npar=4, lhfcalc=True)
        assert calc.int_params["ncore"] == 1
        assert calc.int_params["npar"] is None


def test_ismear():
    atoms = bulk("Cu")

    calc = Vasp(atoms, nsw=10)
    assert calc.int_params["ismear"] is None

    calc = Vasp(atoms, nsw=0)
    assert calc.int_params["ismear"] is None

    calc = Vasp(atoms, kpts=(10, 10, 10), nsw=0)
    assert calc.int_params["ismear"] == -5


def test_ismear_aggressive():
    with change_settings({"VASP_INCAR_COPILOT": "aggressive"}):
        atoms = bulk("Cu")

        calc = Vasp(atoms, nsw=10)
        assert calc.int_params["ismear"] is None

        calc = Vasp(atoms, ismear=-5, nsw=10)
        assert calc.int_params["ismear"] == 1
        assert calc.float_params["sigma"] == 0.1

        calc = Vasp(atoms, ismear=-5, nsw=0)
        assert calc.int_params["ismear"] == 0

        calc = Vasp(atoms, kpts=(10, 10, 10), ismear=-5, nsw=0)
        assert calc.int_params["ismear"] == -5

        calc = Vasp(atoms, ismear=0, nsw=10)
        assert calc.int_params["ismear"] == 0

        calc = Vasp(atoms, nsw=0)
        assert calc.int_params["ismear"] is None

        calc = Vasp(atoms, ismear=-5, nsw=0)
        assert calc.int_params["ismear"] == 0

        calc = Vasp(atoms, kpts=(10, 10, 10), nsw=0)
        assert calc.int_params["ismear"] == -5

        calc = Vasp(atoms, pmg_kpts={"line_density": 100}, ismear=1)
        assert calc.int_params["ismear"] == 0
        assert calc.float_params["sigma"] == 0.01
        assert calc.todict()["reciprocal"] is True

        calc = Vasp(atoms, pmg_kpts={"line_density": 100}, ismear=0, sigma=1e-3)
        assert calc.int_params["ismear"] == 0
        assert calc.float_params["sigma"] == 1e-3
        assert calc.todict()["reciprocal"] is True

        calc = Vasp(atoms, pmg_kpts={"line_density": 100}, ismear=-5)
        assert calc.int_params["ismear"] == 0
        assert calc.todict()["reciprocal"] is True

        calc = Vasp(atoms, kspacing=1.0, ismear=-5)
        assert calc.int_params["ismear"] == 0

        calc = Vasp(atoms, nsw=0, kspacing=1.0, ismear=1, sigma=0.1)
        assert calc.int_params["ismear"] == 1
        assert calc.float_params["sigma"] == 0.1

        atoms[0].symbol = "H"
        calc = Vasp(atoms, kpts=(10, 10, 10), ismear=-5, nsw=10)
        assert calc.int_params["ismear"] == -5


def test_laechg_aggressive():
    with change_settings({"VASP_INCAR_COPILOT": "aggressive"}):
        atoms = bulk("Cu")
        calc = Vasp(atoms, nsw=10, laechg=True)
        assert not calc.bool_params["laechg"]

        calc = Vasp(atoms, laechg=True)
        assert calc.bool_params["laechg"]

        calc = Vasp(atoms, nsw=0, laechg=True)
        assert calc.bool_params["laechg"]

        calc = Vasp(atoms, nsw=0, laechg=False)
        assert not calc.bool_params["laechg"]


def test_ldauprint():
    atoms = bulk("Cu")

    calc = Vasp(atoms, ldau=True)
    assert calc.int_params["ldauprint"] == 1

    calc = Vasp(atoms, ldau_luj={"Cu": {"L": 2, "U": 5, "J": 0.0}})
    assert calc.int_params["ldauprint"] == 1


def test_ldauprint_aggressive():
    with change_settings({"VASP_INCAR_COPILOT": "aggressive"}):
        atoms = bulk("Cu")

        calc = Vasp(atoms, ldau=True)
        assert calc.int_params["ldauprint"] == 1

        calc = Vasp(atoms, ldau=True, ldauprint=0)
        assert calc.int_params["ldauprint"] == 1

        calc = Vasp(atoms, ldau=False, ldauprint=1)
        assert calc.int_params["ldauprint"] is None

        calc = Vasp(atoms, ldau_luj={"Cu": {"L": 2, "U": 5, "J": 0.0}})
        assert calc.int_params["ldauprint"] == 1


def test_lreal():
    atoms = bulk("Cu")

    calc = Vasp(atoms, nsw=10)
    assert calc.special_params["lreal"] is None


def test_lreal_aggressive():
    with change_settings({"VASP_INCAR_COPILOT": "aggressive"}):
        atoms = bulk("Cu")
        calc = Vasp(atoms, lreal=True, nsw=0)
        assert calc.special_params["lreal"] is False

        calc = Vasp(atoms, lreal=True, nsw=10)
        assert calc.special_params["lreal"] is False

        calc = Vasp(atoms, nsw=10)
        assert calc.special_params["lreal"] is None

        calc = Vasp(atoms, lreal="auto", nsw=10)
        assert calc.special_params["lreal"] is False

        calc = Vasp(atoms * (4, 4, 4), lreal="auto", nsw=10)
        assert calc.special_params["lreal"] == "auto"

        calc = Vasp(atoms * (4, 4, 4), lreal=False, nsw=10)
        assert calc.special_params["lreal"] is False


def test_lorbit():
    atoms = bulk("Cu")
    calc = Vasp(atoms, ispin=2)
    assert calc.int_params["lorbit"] == 11

    atoms = bulk("Cu")
    calc = Vasp(atoms, ispin=1)
    assert calc.int_params["lorbit"] is None

    atoms = bulk("Cu")
    atoms.set_initial_magnetic_moments([1.0] * len(atoms))
    calc = Vasp(atoms)
    assert calc.int_params["lorbit"] == 11


def test_dispersion():
    atoms = bulk("Cu")
    calc = Vasp(atoms, xc="r2scan", ivdw=13)
    assert calc.float_params["vdw_s6"] == 1.0
    assert calc.float_params["vdw_s8"] == 0.60187490
    assert calc.float_params["vdw_a1"] == 0.51559235
    assert calc.float_params["vdw_a2"] == 5.77342911

    calc = Vasp(atoms, metagga="r2scan", ivdw=13)
    assert calc.float_params["vdw_s6"] == 1.0
    assert calc.float_params["vdw_s8"] == 0.60187490
    assert calc.float_params["vdw_a1"] == 0.51559235
    assert calc.float_params["vdw_a2"] == 5.77342911

    calc = Vasp(atoms, metagga="r2scan", vdw_s6=2.0)
    assert calc.float_params["vdw_s6"] == 2.0

    calc = Vasp(atoms, metagga="r2scan")
    assert calc.float_params["vdw_s6"] is None

    calc = Vasp(atoms, xc="hse06")
    assert calc.float_params["vdw_s6"] is None

    calc = Vasp(atoms, lhfcalc=True, hfscreen=0.25, ivdw=12)
    assert calc.float_params["vdw_s6"] is None

    calc = Vasp(atoms, xc="hse06", ivdw=12)
    assert calc.float_params["vdw_s6"] == 1.0
    assert calc.float_params["vdw_s8"] == 2.310
    assert calc.float_params["vdw_a1"] == 0.383
    assert calc.float_params["vdw_a2"] == 5.685

    calc = Vasp(atoms, lhfcalc=True, hfscreen=0.2, ivdw=12)
    assert calc.float_params["vdw_s6"] == 1.0
    assert calc.float_params["vdw_s8"] == 2.310
    assert calc.float_params["vdw_a1"] == 0.383
    assert calc.float_params["vdw_a2"] == 5.685


def test_setups():
    atoms = bulk("Cu")
    calc = Vasp(atoms, preset="BulkSet")
    assert calc.parameters["setups"]["Cu"] == ""

    atoms = bulk("Cu")
    calc = Vasp(atoms, preset="SlabSet")
    assert calc.parameters["setups"]["Ba"] == "_sv"
    assert calc.parameters["setups"]["Cu"] == ""

    atoms = bulk("Cu")
    parameters = MPtoASEConverter(atoms=atoms).convert_dict_set(MPScanRelaxSet)
    calc = Vasp(atoms, **parameters)
    assert calc.parameters["setups"]["Cu"] == "_pv"

    atoms = bulk("Cu")
    calc = Vasp(atoms, setups=FILE_DIR / "test_setups.yaml", preset="BulkSet")
    assert calc.parameters["setups"]["Cu"] == "_sv"

    atoms = bulk("Cu")
    calc = Vasp(atoms, setups="setups_recommended.yaml", preset="BulkSet")
    assert calc.parameters["setups"]["Cu"] == ""

    atoms = bulk("Cu")
    calc = Vasp(atoms, setups="minimal", preset="BulkSet")
    assert isinstance(calc.parameters["setups"], str)
    assert calc.parameters["setups"] == "minimal"

    atoms = bulk("Cu")
    calc = Vasp(atoms, preset="QMOFSet")
    assert calc.parameters["setups"]["Cu"] == ""
    assert calc.parameters["setups"]["Er"] == "_3"
    assert calc.parameters["setups"]["Yb"] == "_3"


def test_kpoint_schemes():
    atoms = bulk("Cu")
    calc = Vasp(atoms, kpts=[1, 1, 1], preset="BulkSet")
    assert calc.kpts == [1, 1, 1]

    atoms = bulk("Cu")
    calc = Vasp(atoms, pmg_kpts={"kppa": 1000}, gamma=False)
    assert calc.kpts == [10, 10, 10]
    assert calc.input_params["gamma"] is False

    atoms = bulk("Cu")
    calc = Vasp(atoms, pmg_kpts={"kppa": 1000})
    assert calc.kpts == [10, 10, 10]
    assert calc.input_params["gamma"] is True
    assert isinstance(calc.kpts[0], int)

    atoms = bulk("Cu")
    calc = Vasp(atoms, preset="BulkSet", pmg_kpts={"kppa": 1000}, gamma=False)
    assert calc.kpts == [10, 10, 10]
    assert calc.input_params["gamma"] is False

    atoms = bulk("Cu")
    calc = Vasp(atoms, pmg_kpts={"kppa": 1000}, gamma=True)
    assert calc.kpts == [10, 10, 10]
    assert calc.input_params["gamma"] is True

    atoms = bulk("Cu")
    calc = Vasp(atoms, pmg_kpts={"kppvol": 100})
    assert calc.kpts == [12, 12, 12]

    atoms = bulk("Cu")
    calc = Vasp(atoms, pmg_kpts={"kppvol": 100, "kppa": 1000})
    assert calc.kpts == [12, 12, 12]

    atoms = bulk("Cu")
    calc = Vasp(atoms, pmg_kpts={"kppvol": 10, "kppa": 1000})
    assert calc.kpts == [10, 10, 10]

    atoms = bulk("Cu")
    calc = Vasp(atoms, pmg_kpts={"length_densities": [50, 50, 1]})
    assert calc.kpts == [20, 20, 1]

    atoms = bulk("Cu")
    calc = Vasp(atoms, pmg_kpts={"line_density": 100})
    assert calc.todict()["reciprocal"] is True
    assert calc.kpts[-1] == pytest.approx(np.array([0.375, 0.75, 0.375]))


def test_constraints():
    atoms = bulk("Cu")
    atoms.set_constraint(FixAtoms(indices=[0]))
    calc = Vasp(atoms)
    atoms.calc = calc
    assert isinstance(atoms.constraints[0], FixAtoms)

    atoms = bulk("Cu") * (2, 1, 1)
    atoms.set_constraint(FixBondLength(0, 1))
    with pytest.raises(
        ValueError, match="Atoms object has a constraint that is not compatible"
    ):
        calc = Vasp(atoms)


def test_envvars():
    with change_settings(
        {
            "VASP_PP_PATH": str(Path("/path/to/pseudos")),
            "VASP_VDW": str(Path("/path/to/kernel")),
        }
    ):
        atoms = bulk("Cu")
        atoms.calc = Vasp(atoms, xc="beef-vdw")
        assert os.environ.get("VASP_PP_PATH") == str(Path("/path/to/pseudos"))
        assert os.environ.get("ASE_VASP_VDW") == str(Path("/path/to/kernel"))


def test_bad():
    atoms = bulk("Cu")

    with pytest.raises(ValueError, match="Unsupported k-point generation scheme"):
        Vasp(atoms, pmg_kpts={"test": [100]})

    with pytest.raises(FileNotFoundError):
        Vasp(atoms, preset="BadRelaxSet")


def test_preset_override():
    atoms = bulk("Cu")

    calc = Vasp(atoms, preset="BulkSet", efermi=None)
    assert calc.parameters.get("efermi") is None


def test_logging(caplog):
    atoms = bulk("Cu")
    with caplog.at_level(INFO):
        Vasp(atoms, nsw=0, kpts=(3, 3, 3))
    assert "Recommending LMAXMIX = 4" in caplog.text
    assert "Recommending ISMEAR = -5" in caplog.text
    assert "ismear': -5" in caplog.text
    assert "lmaxmix': 4" in caplog.text

    with caplog.at_level(INFO):
        Vasp(atoms, nsw=0, kpts=(2, 2, 1), ismear=0)
    assert "Recommending LMAXMIX = 4" in caplog.text
    assert "Recommending ISMEAR = -5" in caplog.text
    assert "lmaxmix': 4" in caplog.text
    assert "ismear: -5" not in caplog.text


def test_bad_pmg_converter():
    with pytest.raises(ValueError, match="Either atoms or prev_dir must be provided"):
        MPtoASEConverter()


def test_pmg_input_set():
    atoms = bulk("Cu")
    parameters = MPtoASEConverter(atoms=atoms).convert_dict_set(MPRelaxSet)
    calc = Vasp(atoms, incar_copilot="off", **parameters)
    assert calc.parameters == {
        "algo": "fast",
        "ediff": 5e-05,
        "encut": 520,
        "ibrion": 2,
        "isif": 3,
        "ismear": -5,
        "ispin": 2,
        "lasph": True,
        "lorbit": 11,
        "lreal": "auto",
        "lwave": False,
        "nelm": 100,
        "nsw": 99,
        "pp": "pbe",
        "prec": "accurate",
        "sigma": 0.05,
        "magmom": [0.6],
        "lmaxmix": 4,
        "kpts": (11, 11, 11),
        "gamma": True,
        "setups": {"Cu": "_pv"},
    }


def test_pmg_input_set2():
    atoms = bulk("Fe") * (2, 1, 1)
    atoms[0].symbol = "O"
    parameters = MPtoASEConverter(atoms=atoms).convert_dict_set(MPRelaxSet)
    calc = Vasp(atoms, incar_copilot="off", **parameters)
    assert calc.parameters == {
        "algo": "fast",
        "ediff": 0.0001,
        "encut": 520,
        "ibrion": 2,
        "isif": 3,
        "ismear": -5,
        "ispin": 2,
        "lasph": True,
        "ldau": True,
        "ldauj": [0, 0],
        "ldaul": [0, 2.0],
        "ldautype": 2,
        "ldauu": [0, 5.3],
        "ldauprint": 1,
        "lorbit": 11,
        "lreal": "auto",
        "lwave": False,
        "nelm": 100,
        "nsw": 99,
        "pp": "pbe",
        "prec": "accurate",
        "sigma": 0.05,
        "magmom": [2.3, 2.3],
        "lmaxmix": 4,
        "kpts": (5, 11, 11),
        "gamma": True,
        "setups": {"Fe": "_pv", "O": ""},
    }


def test_ldau_mp():
    atoms = Atoms("HOHMn")
    atoms.center(vacuum=10)
    atoms.pbc = True
    atoms[0].position += 0.1
    atoms[1].position -= 0.1
    parameters = MPtoASEConverter(atoms=atoms).convert_dict_set(MPRelaxSet)
    assert len(parameters["ldauu"]) == 3
    assert parameters["ldauu"] == [0, 0, 3.9]
    assert len(parameters["magmom"]) == 4
    assert parameters["magmom"] == [0.6, 0.6, 0.6, 5.0]


def test_command():
    atoms = bulk("Cu")
    assert Vasp(atoms, kpts=[1, 1, 1]).command.strip() == "vasp_gam"
    assert Vasp(atoms, kpts=[2, 1, 1]).command.strip() == "vasp_std"
    assert Vasp(atoms, kspacing=50).command.strip() == "vasp_gam"
    assert Vasp(atoms, kspacing=0.1).command.strip() == "vasp_std"


@pytest.mark.skipif(which(get_settings().VASP_CMD), reason="VASP is installed")
def test_run(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Cu")
    calc = Vasp(atoms, xc="PBE", use_custodian=False)
    assert calc._run()[0] > 0

    atoms = bulk("Cu")
    calc = Vasp(atoms, xc="PBE", use_custodian=True)
    with pytest.raises(FileNotFoundError):
        calc._run()
