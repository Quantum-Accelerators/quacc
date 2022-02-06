import pytest
from ase.atoms import Atoms
import numpy as np


@pytest.fixture(autouse=True)
def set_env(monkeypatch):
    # Set dummy environment variables in pytest so the tests can run
    monkeypatch.setenv("VASP_PARALLEL_CMD", "")
    monkeypatch.setenv("VASP_PP_PATH", ".")


def mock_get_potential_energy(self):
    # Instead of running .get_potential_energy(), we mock it by attaching
    # dummy results to the atoms object and returning a fake energy. This
    # works because the real atoms.get_potential_energy() takes one argument
    # (i.e. self) and that self object is the atoms object.
    e = -1.0
    self.calc.results = {
        "energy": e,
        "magmom": 0.0,
        "magmoms": [0.0] * len(self),
        "free_energy": -1.0,
        "forces": np.array([[0.0, 0.0, 0.0]] * len(self)),
    }
    return e


@pytest.fixture(autouse=True)
def patch_get_potential_energy(monkeypatch):
    # Monkeypatch the .get_potential_energy() method of the .Atoms object so
    # we aren't running the actual calculation during testing.
    monkeypatch.setattr(Atoms, "get_potential_energy", mock_get_potential_energy)
