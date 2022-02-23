import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms, FixBondLength
from ase.io import read

from quacc.calculators.vasp import SmartVasp
from quacc.defaults.calcs import vasp as v
from quacc.util.atoms import prep_next_run

FILE_DIR = Path(__file__).resolve().parent
DEFAULT_CALCS_DIR = os.path.dirname(v.__file__)
ATOMS_MAG = read(os.path.join(FILE_DIR, "OUTCAR_mag.gz"))
ATOMS_NOMAG = read(os.path.join(FILE_DIR, "OUTCAR_nomag.gz"))
ATOMS_NOSPIN = read(os.path.join(FILE_DIR, "OUTCAR_nospin.gz"))


def test_vanilla_smartvasp():
    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, incar_copilot=False)
    if atoms.calc.asdict() != Vasp().asdict():
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, custodian=False, incar_copilot=False)
    if atoms.calc.asdict() != Vasp().asdict():
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, encut=None, incar_copilot=False)
    if atoms.calc.asdict() != Vasp().asdict():
        raise AssertionError


def test_presets():
    atoms = bulk("Co") * (2, 2, 1)
    atoms[-1].symbol = "Fe"

    atoms_fullpath = SmartVasp(
        atoms, preset=os.path.join(DEFAULT_CALCS_DIR, "BulkRelaxSet")
    )
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms != atoms_fullpath:
        raise AssertionError
    if atoms.calc.xc.lower() != "pbe":
        raise AssertionError
    if atoms.calc.string_params["algo"] != "fast":
        raise AssertionError
    if atoms.calc.exp_params["ediff"] != 1e-5:
        raise AssertionError
    if atoms.calc.int_params["isif"] != 3:
        raise AssertionError

    atoms = SmartVasp(atoms, xc="rpbe", preset="SlabRelaxSet")
    if atoms.calc.xc.lower() != "rpbe":
        raise AssertionError
    if atoms.calc.string_params["algo"] != "fast":
        raise AssertionError
    if atoms.calc.exp_params["ediff"] != 1e-5:
        raise AssertionError
    if atoms.calc.int_params["isif"] != 2:
        raise AssertionError

    atoms = SmartVasp(atoms, xc="scan", preset="MPScanRelaxSet")
    if atoms.calc.xc.lower() != "scan":
        raise AssertionError
    if atoms.calc.string_params["algo"] != "all":
        raise AssertionError
    if atoms.calc.exp_params["ediff"] != 1e-5:
        raise AssertionError
    if atoms.calc.int_params["isif"] != 3:
        raise AssertionError


def test_lmaxmix():
    atoms = bulk("Cu")
    atoms = SmartVasp(atoms)
    if atoms.calc.int_params["lmaxmix"] != 4:
        raise AssertionError

    atoms = bulk("Ce")
    atoms = SmartVasp(atoms)
    if atoms.calc.int_params["lmaxmix"] != 6:
        raise AssertionError

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[-1].symbol = "Ce"
    atoms = SmartVasp(atoms)
    if atoms.calc.int_params["lmaxmix"] != 6:
        raise AssertionError


def test_autodipole():
    atoms = bulk("Cu")
    com = atoms.get_center_of_mass(scaled=True)
    atoms = SmartVasp(atoms, auto_dipole=True)
    if atoms.calc.bool_params["ldipol"] is not True:
        raise AssertionError
    if atoms.calc.int_params["idipol"] != 3:
        raise AssertionError
    if not np.array_equal(atoms.calc.list_float_params["dipol"], com):
        raise AssertionError

    atoms = SmartVasp(atoms, auto_dipole=True, idipol=2)
    if atoms.calc.bool_params["ldipol"] is not True:
        raise AssertionError
    if atoms.calc.int_params["idipol"] != 2:
        raise AssertionError
    if not np.array_equal(atoms.calc.list_float_params["dipol"], com):
        raise AssertionError

    atoms = SmartVasp(atoms, preset="SlabRelaxSet")
    if atoms.calc.bool_params["ldipol"] is not True:
        raise AssertionError
    if atoms.calc.int_params["idipol"] != 3:
        raise AssertionError
    if not np.array_equal(atoms.calc.list_float_params["dipol"], com):
        raise AssertionError

    atoms = SmartVasp(atoms, preset="SlabRelaxSet", idipol=2)
    if atoms.calc.bool_params["ldipol"] is not True:
        raise AssertionError
    if atoms.calc.int_params["idipol"] != 2:
        raise AssertionError
    if not np.array_equal(atoms.calc.list_float_params["dipol"], com):
        raise AssertionError

    atoms = SmartVasp(atoms, auto_dipole=False, preset="SlabRelaxSet", idipol=2)
    if atoms.calc.bool_params["ldipol"] is not None:
        raise AssertionError
    if atoms.calc.int_params["idipol"] != 2:
        raise AssertionError
    if atoms.calc.list_float_params["dipol"] is not None:
        raise AssertionError


def test_kspacing():
    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, kspacing=0.1, ismear=-5)
    if atoms.calc.int_params["ismear"] != -5:
        raise AssertionError

    atoms = SmartVasp(atoms, kspacing=100, ismear=-5)
    if atoms.calc.int_params["ismear"] != 0:
        raise AssertionError


def test_magmoms():
    atoms = bulk("Mg")
    atoms = SmartVasp(atoms)
    if atoms.has("initial_magmoms") is not False:
        raise AssertionError

    atoms = bulk("Mg")
    atoms.set_initial_magnetic_moments([3.14] * len(atoms))
    atoms = SmartVasp(atoms)
    if atoms.get_initial_magnetic_moments().tolist() != [3.14] * len(atoms):
        raise AssertionError

    atoms = bulk("Cu") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.get_initial_magnetic_moments().tolist() != [2.0] * (len(atoms) - 1) + [
        5.0
    ]:
        raise AssertionError

    atoms = bulk("Zn") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    atoms = SmartVasp(atoms, preset="BulkRelaxSet", mag_default=2.5)
    if atoms.get_initial_magnetic_moments().tolist() != [2.5] * (len(atoms) - 1) + [
        5.0
    ]:
        raise AssertionError

    atoms = bulk("Eu") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    atoms = SmartVasp(atoms, preset="SlabRelaxSet")
    if atoms.get_initial_magnetic_moments().tolist() != [7.0] * (len(atoms) - 1) + [
        5.0
    ]:
        raise AssertionError

    atoms = bulk("Cu") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    atoms = SmartVasp(atoms, preset="MPScanRelaxSet")
    if atoms.get_initial_magnetic_moments().tolist() != [1.0] * (len(atoms) - 1) + [
        5.0
    ]:
        raise AssertionError

    atoms = bulk("Cu") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    atoms.set_initial_magnetic_moments([3.14] * (len(atoms) - 1) + [1.0])
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.get_initial_magnetic_moments().tolist() != [3.14] * (
        len(atoms) - 1
    ) + [1.0]:
        raise AssertionError

    atoms = bulk("Co") * (2, 2, 1)
    atoms[-1].symbol = "Fe"
    atoms.set_initial_magnetic_moments([0.0] * len(atoms))
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.has("initial_magmoms") is not True:
        raise AssertionError
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError

    atoms = deepcopy(ATOMS_MAG)
    mags = atoms.get_magnetic_moments()
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if np.array_equal(atoms.get_initial_magnetic_moments(), mags) is not True:
        raise AssertionError

    atoms = deepcopy(ATOMS_MAG)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet", mag_cutoff=2.0)
    if atoms.has("initial_magmoms") is not True:
        raise AssertionError
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError

    atoms = deepcopy(ATOMS_MAG)
    if atoms.get_magnetic_moments()[0] != 0.468:
        raise AssertionError
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.get_initial_magnetic_moments()[0] != 0.468:
        raise AssertionError

    atoms = deepcopy(ATOMS_NOMAG)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.has("initial_magmoms") is not True:
        raise AssertionError
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError
    atoms.calc.results = {"energy": -1.0}  # mock calculation run
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.has("initial_magmoms") is not True:
        raise AssertionError
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError

    atoms = deepcopy(ATOMS_NOMAG)
    calc = SinglePointDFTCalculator(atoms, **{"magmoms": [0.0] * len(atoms)})
    atoms.calc = calc
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError

    atoms = deepcopy(ATOMS_NOMAG)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError
    atoms.calc.results = {"magmoms": [0.0] * len(atoms)}
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError

    atoms = deepcopy(ATOMS_NOMAG)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.has("initial_magmoms") is not True:
        raise AssertionError
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError
    atoms.calc.results = {"magmoms": [1.0] * len(atoms)}
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == 1):
        raise AssertionError

    atoms = deepcopy(ATOMS_NOSPIN)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.has("initial_magmoms") is not True:
        raise AssertionError
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError
    atoms.calc.results = {"magmoms": [1.0] * len(atoms)}
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == 1):
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, preset="BulkRelaxSet", copy_magmoms=False)
    if not np.all(atoms.get_initial_magnetic_moments() == 2.0):
        raise AssertionError
    atoms.calc.results = {"magmoms": [3.0] * len(atoms)}
    atoms = SmartVasp(atoms, preset="BulkRelaxSet", copy_magmoms=False)
    if not np.all(atoms.get_initial_magnetic_moments() == 2.0):
        raise AssertionError

    atoms = bulk("Mg")
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == 1.0):
        raise AssertionError
    atoms.calc.results = {"magmoms": [0.0] * len(atoms)}
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError

    atoms = bulk("Mg")
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == 1.0):
        raise AssertionError
    atoms.calc.results = {"magmoms": [-0.01] * len(atoms)}
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError

    atoms = bulk("Mg")
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == 1.0):
        raise AssertionError
    atoms.calc.results = {"magmoms": [-5] * len(atoms)}
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if not np.all(atoms.get_initial_magnetic_moments() == -5):
        raise AssertionError

    atoms = deepcopy(ATOMS_MAG)
    mags = atoms.get_magnetic_moments()
    atoms = prep_next_run(atoms)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.has("initial_magmoms") is not True:
        raise AssertionError
    if atoms.get_initial_magnetic_moments().tolist() != mags.tolist():
        raise AssertionError

    atoms = deepcopy(ATOMS_NOMAG)
    atoms = prep_next_run(atoms)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.has("initial_magmoms") is not True:
        raise AssertionError
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError

    atoms = deepcopy(ATOMS_NOSPIN)
    atoms = prep_next_run(atoms)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.has("initial_magmoms") is not True:
        raise AssertionError
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError

    atoms = deepcopy(ATOMS_MAG)
    atoms = prep_next_run(atoms)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet", mag_cutoff=10.0)
    if atoms.has("initial_magmoms") is not True:
        raise AssertionError
    if not np.all(atoms.get_initial_magnetic_moments() == 0):
        raise AssertionError

    atoms = bulk("Mg")
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"energy": -1.0, "magmoms": [0.0] * len(atoms)}
    atoms = prep_next_run(atoms)
    atoms *= (2, 2, 2)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    np.all(atoms.get_initial_magnetic_moments() == 0.0)

    atoms = bulk("Mg")
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"energy": -1.0, "magmoms": [-0.02] * len(atoms)}
    atoms = prep_next_run(atoms)
    atoms *= (2, 2, 2)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    np.all(atoms.get_initial_magnetic_moments() == 0.0)

    atoms = bulk("Mg")
    atoms = SmartVasp(atoms)
    atoms.calc.results = {"energy": -1.0}
    atoms = prep_next_run(atoms)
    atoms *= (2, 2, 2)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    np.all(atoms.get_initial_magnetic_moments() == 0.0)

    atoms = bulk("Mg")
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    atoms.calc.results = {"energy": -1.0, "magmoms": [3.14] * len(atoms)}
    atoms = prep_next_run(atoms)
    atoms *= (2, 2, 2)
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    np.all(atoms.get_initial_magnetic_moments() == 3.14)


def test_unused_flags():
    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, preset="BulkRelaxSet", potim=1.5, nsw=0)
    if atoms.calc.int_params["nsw"] != 0:
        raise AssertionError
    if atoms.calc.exp_params["ediffg"] is not None:
        raise AssertionError
    if atoms.calc.int_params["isif"] is not None:
        raise AssertionError
    if atoms.calc.float_params["potim"] is not None:
        raise AssertionError

    atoms = SmartVasp(atoms, ldau=False, ldauprint=2)
    if atoms.calc.int_params["ldauprint"] is not None:
        raise AssertionError
    if atoms.calc.bool_params["ldau"] is not None:
        raise AssertionError


def test_lasph():
    atoms = bulk("Cu")

    atoms = SmartVasp(atoms, xc="rpbe")
    if atoms.calc.bool_params["lasph"] is not None:
        raise AssertionError

    atoms = SmartVasp(atoms, xc="m06l")
    if atoms.calc.bool_params["lasph"] is not True:
        raise AssertionError

    atoms = SmartVasp(atoms, xc="m06l", lasph=False)
    if atoms.calc.bool_params["lasph"] is not True:
        raise AssertionError

    atoms = SmartVasp(atoms, xc="hse06")
    if atoms.calc.bool_params["lasph"] is not True:
        raise AssertionError

    atoms = SmartVasp(atoms, xc="beef-vdw")
    if atoms.calc.bool_params["lasph"] is not True:
        raise AssertionError

    atoms = SmartVasp(atoms, ldau_luj={"Cu": {"L": 2, "U": 5, "J": 0.0}})
    if atoms.calc.bool_params["lasph"] is not True:
        raise AssertionError


def test_lmaxtau():
    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, lasph=True)
    if atoms.calc.int_params["lmaxtau"] is not None:
        raise AssertionError

    atoms = bulk("Ce")
    atoms = SmartVasp(atoms, lasph=True)
    if atoms.calc.int_params["lmaxtau"] != 8:
        raise AssertionError

    atoms = bulk("Cu") * (2, 2, 2)
    atoms[-1].symbol = "Ce"
    atoms = SmartVasp(atoms, lasph=True)
    if atoms.calc.int_params["lmaxtau"] != 8:
        raise AssertionError


def test_algo():
    atoms = bulk("Cu")

    atoms = SmartVasp(atoms, xc="rpbe")
    if atoms.calc.string_params["algo"] is not None:
        raise AssertionError

    atoms = SmartVasp(atoms, xc="m06l")
    if atoms.calc.string_params["algo"] != "all":
        raise AssertionError

    atoms = SmartVasp(atoms, xc="m06l", algo="fast")
    if atoms.calc.string_params["algo"] != "all":
        raise AssertionError

    atoms = SmartVasp(atoms, xc="hse06")
    if atoms.calc.string_params["algo"] != "damped":
        raise AssertionError
    if atoms.calc.float_params["time"] != 0.5:
        raise AssertionError

    atoms[0].symbol = "H"
    atoms = SmartVasp(atoms, xc="hse06")
    if atoms.calc.string_params["algo"] != "all":
        raise AssertionError


def test_kpar():
    atoms = bulk("Cu")

    atoms = SmartVasp(atoms, kpts=[2, 2, 1], kpar=4)
    if atoms.calc.int_params["kpar"] != 4:
        raise AssertionError

    atoms = SmartVasp(atoms, kpar=4)
    if atoms.calc.int_params["kpar"] != 1:
        raise AssertionError


def test_isym():
    atoms = bulk("Cu")

    atoms = SmartVasp(atoms, isym=2)
    if atoms.calc.int_params["isym"] != 2:
        raise AssertionError

    atoms = SmartVasp(atoms, isym=0)
    if atoms.calc.int_params["isym"] != 0:
        raise AssertionError

    atoms = SmartVasp(atoms, xc="hse06", isym=2)
    if atoms.calc.int_params["isym"] != 3:
        raise AssertionError

    atoms = SmartVasp(atoms, isym=2, nsw=100)
    if atoms.calc.int_params["isym"] != 0:
        raise AssertionError

    atoms = SmartVasp(atoms, xc="hse06", isym=2, nsw=100)
    if atoms.calc.int_params["isym"] != 0:
        raise AssertionError


def test_ncore():
    atoms = bulk("Cu")

    atoms = SmartVasp(atoms, ncore=16)
    if atoms.calc.int_params["ncore"] != 1:
        raise AssertionError

    atoms = SmartVasp(atoms, npar=16)
    if atoms.calc.int_params["ncore"] != 1:
        raise AssertionError
    if atoms.calc.int_params["npar"] is not None:
        raise AssertionError

    atoms *= (2, 2, 2)
    atoms = SmartVasp(atoms, ncore=4)
    if atoms.calc.int_params["ncore"] != 4:
        raise AssertionError

    atoms = SmartVasp(atoms, ncore=4, lhfcalc=True)
    if atoms.calc.int_params["ncore"] != 1:
        raise AssertionError

    atoms = SmartVasp(atoms, npar=4, lhfcalc=True)
    if atoms.calc.int_params["ncore"] != 1:
        raise AssertionError
    if atoms.calc.int_params["npar"] is not None:
        raise AssertionError


def test_ismear():
    atoms = bulk("Cu")

    atoms = SmartVasp(atoms, nsw=10)
    if atoms.calc.int_params["ismear"] is not None:
        raise AssertionError

    atoms = SmartVasp(atoms, ismear=-5, nsw=10)
    if atoms.calc.int_params["ismear"] != 1:
        raise AssertionError
    if atoms.calc.float_params["sigma"] != 0.1:
        raise AssertionError

    atoms = SmartVasp(atoms, ismear=-5, nsw=0)
    if atoms.calc.int_params["ismear"] != 0:
        raise AssertionError

    atoms = SmartVasp(atoms, kpts=(10, 10, 10), ismear=-5, nsw=0)
    if atoms.calc.int_params["ismear"] != -5:
        raise AssertionError

    atoms = SmartVasp(atoms, ismear=0, nsw=10)
    if atoms.calc.int_params["ismear"] != 0:
        raise AssertionError

    atoms = SmartVasp(atoms, nedos=3001, nsw=0)
    if atoms.calc.int_params["ismear"] != 0:
        raise AssertionError

    atoms = SmartVasp(atoms, ismear=-5, nedos=3001, nsw=0)
    if atoms.calc.int_params["ismear"] != 0:
        raise AssertionError

    atoms = SmartVasp(atoms, kpts=(10, 10, 10), nedos=3001, nsw=0)
    if atoms.calc.int_params["ismear"] != -5:
        raise AssertionError

    atoms = SmartVasp(atoms, auto_kpts={"line_density": 100}, ismear=1)
    if atoms.calc.int_params["ismear"] != 0:
        raise AssertionError
    if atoms.calc.float_params["sigma"] != 0.01:
        raise AssertionError

    atoms = SmartVasp(atoms, auto_kpts={"line_density": 100}, ismear=0, sigma=1e-3)
    if atoms.calc.int_params["ismear"] != 0:
        raise AssertionError
    if atoms.calc.float_params["sigma"] != 1e-3:
        raise AssertionError

    atoms = SmartVasp(atoms, auto_kpts={"line_density": 100}, ismear=-5)
    if atoms.calc.int_params["ismear"] != 0:
        raise AssertionError

    atoms = SmartVasp(atoms, kspacing=1.0, ismear=-5)
    if atoms.calc.int_params["ismear"] != 0:
        raise AssertionError
    if atoms.calc.float_params["sigma"] != 0.05:
        raise AssertionError

    atoms = SmartVasp(atoms, nsw=0, kspacing=1.0, ismear=1, sigma=0.1)
    if atoms.calc.int_params["ismear"] != 1:
        raise AssertionError
    if atoms.calc.float_params["sigma"] != 0.1:
        raise AssertionError

    atoms[0].symbol = "H"
    atoms = SmartVasp(atoms, kpts=(10, 10, 10), ismear=-5, nsw=10)
    if atoms.calc.int_params["ismear"] != -5:
        raise AssertionError


def test_laechg():
    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, nsw=10, laechg=True)
    if atoms.calc.bool_params["laechg"] is not None:
        raise AssertionError

    atoms = SmartVasp(atoms, laechg=True)
    if atoms.calc.bool_params["laechg"] is not True:
        raise AssertionError

    atoms = SmartVasp(atoms, nsw=0, laechg=True)
    if atoms.calc.bool_params["laechg"] is not True:
        raise AssertionError

    atoms = SmartVasp(atoms, nsw=0, laechg=False)
    if atoms.calc.bool_params["laechg"] is not False:
        raise AssertionError


def test_ldauprint():
    atoms = bulk("Cu")

    atoms = SmartVasp(atoms, ldau=True)
    if atoms.calc.int_params["ldauprint"] != 1:
        raise AssertionError

    atoms = SmartVasp(atoms, ldau=True, ldauprint=0)
    if atoms.calc.int_params["ldauprint"] != 1:
        raise AssertionError

    atoms = SmartVasp(atoms, ldau=False, ldauprint=1)
    if atoms.calc.int_params["ldauprint"] is not None:
        raise AssertionError

    atoms = SmartVasp(atoms, ldau_luj={"Cu": {"L": 2, "U": 5, "J": 0.0}})
    if atoms.calc.int_params["ldauprint"] != 1:
        raise AssertionError


def test_lreal():
    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, lreal=True, nsw=0)
    if atoms.calc.special_params["lreal"] is not False:
        raise AssertionError

    atoms = SmartVasp(atoms, lreal=True, nsw=10)
    if atoms.calc.special_params["lreal"] is not True:
        raise AssertionError

    atoms = SmartVasp(atoms, nsw=10)
    if atoms.calc.special_params["lreal"] is not None:
        raise AssertionError


def test_lorbit():
    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, ispin=2)
    if atoms.calc.int_params["lorbit"] != 11:
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, ispin=1)
    if atoms.calc.int_params["lorbit"] is not None:
        raise AssertionError

    atoms = bulk("Cu")
    atoms.set_initial_magnetic_moments([1.0] * len(atoms))
    atoms = SmartVasp(atoms)
    if atoms.calc.int_params["lorbit"] != 11:
        raise AssertionError


def test_setups():
    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, preset="BulkRelaxSet")
    if atoms.calc.parameters["setups"]["Cu"] != "":
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, preset="SlabRelaxSet")
    if atoms.calc.parameters["setups"]["Ba"] != "_sv":
        raise AssertionError
    if atoms.calc.parameters["setups"]["Cu"] != "":
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, preset="MPScanRelaxSet")
    if atoms.calc.parameters["setups"]["Cu"] != "_pv":
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(
        atoms,
        setups=os.path.join(FILE_DIR, "test_setups.yaml"),
        preset="BulkRelaxSet",
    )
    if atoms.calc.parameters["setups"]["Cu"] != "_pv":
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(
        atoms,
        setups="pbe54_MP.yaml",
        preset="BulkRelaxSet",
    )
    if atoms.calc.parameters["setups"]["Cu"] != "_pv":
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, setups="minimal", preset="MPScanRelaxSet")
    if not (
        isinstance(atoms.calc.parameters["setups"], str)
        and atoms.calc.parameters["setups"] == "minimal"
    ):
        raise AssertionError


def test_kpoint_schemes():
    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, kpts=[1, 1, 1], preset="BulkRelaxSet")
    if atoms.calc.kpts != [1, 1, 1]:
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, auto_kpts={"grid_density": 1000}, gamma=False)
    if atoms.calc.kpts != [10, 10, 10]:
        raise AssertionError
    if atoms.calc.input_params["gamma"] is not False:
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, auto_kpts={"grid_density": 1000})
    if atoms.calc.kpts != [10, 10, 10]:
        raise AssertionError
    if atoms.calc.input_params["gamma"] is not True:
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(
        atoms,
        preset="BulkRelaxSet",
        auto_kpts={"grid_density": 1000},
        gamma=False,
    )
    if atoms.calc.kpts != [10, 10, 10]:
        raise AssertionError
    if atoms.calc.input_params["gamma"] is not False:
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, auto_kpts={"grid_density": 1000}, gamma=True)
    if atoms.calc.kpts != [10, 10, 10]:
        raise AssertionError
    if atoms.calc.input_params["gamma"] is not True:
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, auto_kpts={"reciprocal_density": 100})
    if atoms.calc.kpts != [12, 12, 12]:
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, auto_kpts={"max_mixed_density": [100, 1000]})
    if atoms.calc.kpts != [12, 12, 12]:
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, auto_kpts={"max_mixed_density": [10, 1000]})
    if atoms.calc.kpts != [10, 10, 10]:
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, auto_kpts={"length_density": [50, 50, 1]})
    if atoms.calc.kpts != [20, 20, 1]:
        raise AssertionError

    atoms = bulk("Cu")
    atoms = SmartVasp(atoms, auto_kpts={"line_density": 100})
    if atoms.calc.kpts[-1, :] != pytest.approx(
        np.array([1.30537091e00, 1.11022302e-16, 1.30537091e00])
    ):
        raise AssertionError


def test_constraints():
    atoms = bulk("Cu")
    atoms.set_constraint(FixAtoms(indices=[0]))
    atoms = SmartVasp(atoms)
    if not isinstance(atoms.constraints[0], FixAtoms):
        raise AssertionError

    atoms = bulk("Cu") * (2, 1, 1)
    atoms.set_constraint(FixBondLength(0, 1))
    with pytest.raises(ValueError):
        atoms = SmartVasp(atoms)


def test_bad():
    atoms = bulk("Cu")
    with pytest.raises(ValueError):
        atoms = SmartVasp(atoms, auto_kpts={"max_mixed_density": [100]})

    with pytest.raises(ValueError):
        atoms = SmartVasp(atoms, auto_kpts={"length_density": [100]})

    with pytest.raises(ValueError):
        atoms = SmartVasp(atoms, auto_kpts={"test": [100]})

    with pytest.warns(Warning):
        atoms = SmartVasp(atoms, auto_kpts={"max_mixed_density": [1000, 100]})

    with pytest.raises(ValueError):
        atoms = SmartVasp(atoms, preset="BadRelaxSet")


def test_bad_custodian(monkeypatch):
    monkeypatch.setenv("VASP_CUSTODIAN_SETTINGS", ".")
    atoms = bulk("Cu")
    with pytest.raises(FileNotFoundError):
        atoms = SmartVasp(atoms)
