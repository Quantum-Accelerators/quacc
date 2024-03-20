from __future__ import annotations

import os
from pathlib import Path
from unittest import mock

import pytest

FILE_DIR = Path(__file__).parent
GULP_DIR = Path(FILE_DIR, "gulp_run")


@pytest.fixture(autouse=True)
def mock_settings_env_vars():
    with mock.patch.dict(os.environ, {"GULP_LIB": "."}):
        yield


def mock_execute(self, *args, **kwargs):
    from shutil import copy

    copy(GULP_DIR / "gulp.got.gz", Path(self.directory, "gulp.got.gz"))
    copy(GULP_DIR / "gulp.cif.gz", Path(self.directory, "gulp.cif.gz"))
    copy(GULP_DIR / "gulp.xyz.gz", Path(self.directory, "gulp.xyz.gz"))


@pytest.fixture(autouse=True)
def patch_execute(monkeypatch):
    from ase.calculators.calculator import FileIOCalculator

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute)


def mock_read_results(self, *args, **kwargs):
    from ase.calculators.lj import LennardJones
    from ase.io import read

    atoms = read(Path(self.directory, "gulp.xyz.gz"))
    atoms.calc = LennardJones()
    atoms.get_potential_energy()
    self.calc = atoms.calc
    self.results = atoms.calc.results


@pytest.fixture(autouse=True)
def patch_read_results(monkeypatch):
    from ase.calculators.gulp import GULP

    monkeypatch.setattr(GULP, "read_results", mock_read_results)


def mock_get_opt_state(self, **kwargs):
    # Instead of running GULP and getting the opt state, we'll just mock
    # it as true

    self.get_opt_state = True

    return True


@pytest.fixture(autouse=True)
def patch_get_opt_state(monkeypatch):
    from ase.calculators.gulp import GULP

    # Monkeypatch the .get_opt_state() method of the GULP calculator object so
    # we aren't fetching the actual state
    monkeypatch.setattr(GULP, "get_opt_state", mock_get_opt_state)
