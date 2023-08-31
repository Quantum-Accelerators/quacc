import os
from pathlib import Path
from shutil import copy
from unittest import mock

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.gulp import GULP

FILE_DIR = Path(__file__).resolve().parent
GULP_DIR = os.path.join(FILE_DIR, "gulp_run")


@pytest.fixture(autouse=True)
def mock_settings_env_vars():
    with mock.patch.dict(os.environ, {"GULP_LIB": "."}):
        yield


def mock_get_potential_energy(self, **kwargs):
    # Instead of running .get_potential_energy(), we mock it by attaching
    # dummy results to the atoms object and returning a fake energy. This
    # works because the real atoms.get_potential_energy() takes one argument
    # (i.e. self) and that self object is the atoms object.
    e = -1.0
    self.calc.results = {
        "energy": e,
        "forces": np.array([[0.0, 0.0, 0.0]] * len(self)),
    }
    for f in os.listdir(GULP_DIR):
        copy(os.path.join(GULP_DIR, f), f)
    return e


@pytest.fixture(autouse=True)
def patch_get_potential_energy(monkeypatch):
    # Monkeypatch the .get_potential_energy() method of the Atoms object so
    # we aren't running the actual calculation during testing.
    monkeypatch.setattr(Atoms, "get_potential_energy", mock_get_potential_energy)


def mock_get_opt_state(self, **kwargs):
    # Instead of running GULP and getting the opt state, we'll just mock
    # it as true

    self.get_opt_state = True

    return True


@pytest.fixture(autouse=True)
def patch_get_opt_state(monkeypatch):
    # Monkeypatch the .get_opt_state() method of the GULP calculator object so
    # we aren't fetching the actual state
    monkeypatch.setattr(GULP, "get_opt_state", mock_get_opt_state)
