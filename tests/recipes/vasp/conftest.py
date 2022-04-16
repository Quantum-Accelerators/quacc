import numpy as np
import pytest
from ase.atoms import Atoms
from ase.optimize.optimize import Optimizer

from quacc.schemas.calc import summarize_run as calc_summarize_run


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
    return e


@pytest.fixture(autouse=True)
def patch_get_potential_energy(monkeypatch):
    # Monkeypatch the .get_potential_energy() method of the .Atoms object so
    # we aren't running the actual calculation during testing.
    monkeypatch.setattr(Atoms, "get_potential_energy", mock_get_potential_energy)


def mock_dynrun(*args, **kwargs):
    # Instead of running optimizer.run(), just proceed
    return True


@pytest.fixture(autouse=True)
def patch_dynrun(monkeypatch):
    # Monkeypatch the optimizer.run() method of the .Atoms object so
    # we aren't running the actual calculation during testing.
    monkeypatch.setattr(Optimizer, "run", mock_dynrun)


def mock_summarize_run(atoms, **kwargs):
    # Instead of running the VASP-specific summarize_run(), we mock it with the
    # general calculator schema which does not require VASP files to be
    # in the working directory and will work with pytest.

    prep_next_run = kwargs.get("prep_next_run", True)
    additional_fields = kwargs.get("additional_fields")
    output = calc_summarize_run(
        atoms, prep_next_run=prep_next_run, additional_fields=additional_fields
    )
    output["output"] = {"energy": -1.0}
    return output


@pytest.fixture(autouse=True)
def patch_summarize_run(monkeypatch):
    # Monkeypatch the summarize_run() function so that we aren't relying on real
    # VASP files to be in the working directory during the test.
    monkeypatch.setattr("quacc.recipes.vasp.core.summarize_run", mock_summarize_run)
    monkeypatch.setattr("quacc.recipes.vasp.qmof.summarize_run", mock_summarize_run)
    monkeypatch.setattr("quacc.recipes.vasp.slabs.summarize_run", mock_summarize_run)
