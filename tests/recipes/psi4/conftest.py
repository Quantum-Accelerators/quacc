import pytest


def mock_get_potential_energy(self, **kwargs):
    # Instead of running .get_potential_energy(), we mock it by attaching
    # dummy results to the atoms object and returning a fake energy. This
    # works because the real atoms.get_potential_energy() takes one argument
    # (i.e. self) and that self object is the atoms object.
    self.calc.results["energy"] = -1.0
    return self.calc.results["energy"


@pytest.fixture(autouse=True)
def patch_get_potential_energy(monkeypatch):
    # Monkeypatch the .get_potential_energy() method of the .Atoms object so
    # we aren't running the actual calculation during testing.
    monkeypatch.setattr(Atoms, "get_potential_energy", mock_get_potential_energy)
