from pathlib import Path
from shutil import copy

import pytest
from ase.calculators.calculator import FileIOCalculator
from ase.calculators.gaussian import Gaussian
from ase.calculators.lj import LennardJones
from ase.io import read

FILE_DIR = Path(__file__).parent
GAUSSIAN_DIR = Path(FILE_DIR, "gaussian_run")


def mock_execute(self, *args, **kwargs):
    copy(GAUSSIAN_DIR / "Gaussian.log.gz", "Gaussian.log.gz")


@pytest.fixture(autouse=True)
def patch_execute(monkeypatch):
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute)


def mock_read_results(self, *args, **kwargs):
    atoms = read("Gaussian.com")
    atoms.calc = LennardJones()
    atoms.get_potential_energy()
    self.calc = atoms.calc
    self.results = atoms.calc.results


@pytest.fixture(autouse=True)
def patch_read_results(monkeypatch):
    monkeypatch.setattr(Gaussian, "read_results", mock_read_results)
