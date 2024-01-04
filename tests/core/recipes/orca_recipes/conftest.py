import os
from pathlib import Path
from shutil import copy, which

import numpy as np
import pytest
from ase import Atoms

FILE_DIR = Path(__file__).parent
ORCA_DIR = os.path.join(FILE_DIR, "orca_run")


def mock_get_potential_energy(self, **kwargs):
    # Instead of running .get_potential_energy(), we mock it by attaching
    # dummy results to the atoms object and returning a fake energy. This
    # works because the real atoms.get_potential_energy() takes one argument
    # (i.e. self) and that self object is the atoms object.
    e = -1.0
    self.calc.results = {"energy": e, "forces": np.array([[0.0, 0.0, 0.0]] * len(self))}
    for f in os.listdir(ORCA_DIR):
        copy(os.path.join(ORCA_DIR, f), f)
    return e


@pytest.fixture(autouse=True)
def patch_get_potential_energy(monkeypatch):
    # Monkeypatch the .get_potential_energy() method of the Atoms object so
    # we aren't running the actual calculation during testing.
    from quacc import SETTINGS

    orca_path = which(str(SETTINGS.ORCA_CMD))

    has_orca = bool(orca_path and os.path.getsize(orca_path) > 1024 * 1024)
    if not has_orca:
        monkeypatch.setattr(Atoms, "get_potential_energy", mock_get_potential_energy)
