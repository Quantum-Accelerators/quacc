import pytest
from ase.calculators.emt import EMT
from ase.calculators.onetep import OnetepTemplate
from ase.io import read


def mock_execute(self, *args, **kwargs):
    pass


@pytest.fixture(autouse=True)
def patch_execute(monkeypatch):
    monkeypatch.setattr(OnetepTemplate, "execute", mock_execute)


def mock_read_results(self, directory, *args, **kwargs):
    atoms = read(directory / "onetep.dat")
    atoms.calc = EMT()
    atoms.get_potential_energy()
    return atoms.calc.results


@pytest.fixture(autouse=True)
def patch_read_results(monkeypatch):
    monkeypatch.setattr(OnetepTemplate, "read_results", mock_read_results)
