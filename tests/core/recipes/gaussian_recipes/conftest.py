from pathlib import Path

import pytest

FILE_DIR = Path(__file__).parent
GAUSSIAN_DIR = Path(FILE_DIR, "gaussian_run")


def mock_execute(self, *args, **kwargs):
    from shutil import copy
    copy(GAUSSIAN_DIR / "Gaussian.log.gz", Path(self.directory, "Gaussian.log.gz"))


@pytest.fixture(autouse=True)
def patch_execute(monkeypatch):
    from ase.calculators.calculator import FileIOCalculator

    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute)


def mock_read_results(self, *args, **kwargs):
    from ase.calculators.lj import LennardJones
    from ase.io import read
    atoms = read(Path(self.directory, "Gaussian.com"))
    atoms.calc = LennardJones()
    atoms.get_potential_energy()
    self.calc = atoms.calc
    self.results = atoms.calc.results


@pytest.fixture(autouse=True)
def patch_read_results(monkeypatch):
    from ase.calculators.gaussian import Gaussian

    monkeypatch.setattr(Gaussian, "read_results", mock_read_results)
