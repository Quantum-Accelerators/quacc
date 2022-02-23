import numpy as np
import pytest
from ase.atoms import Atoms

from quacc.schemas.calc import summarize_run as calc_summarize_run


def mock_get_potential_energy(self, **kwargs):
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


def mock_summarize_run(atoms, **kwargs):
    # Instead of running the VASP-specific summarize_run(), we mock it with the
    # general calculator schema which does not require VASP files to be
    # in the working directory and will work with pytest.

    prep_next_run = kwargs.get("prep_next_run", True)
    additional_fields = kwargs.get("additional_fields",None)
    output = calc_summarize_run(
        atoms, prep_next_run=prep_next_run, additional_fields=additional_fields
    )
    return output


@pytest.fixture(autouse=True)
def patch_summarize_run(monkeypatch):
    # Monkeypatch the summarize_run() function so that we aren't relying on real
    # VASP files to be in the working directory during the test. Note that even though
    # summarize_run() is a function in the quacc.schemas.vasp module, we modify it
    # only in quacc.recipes.vasp.core/.slabs because otherwise it will not work properly.
    monkeypatch.setattr("quacc.recipes.vasp.core.summarize_run", mock_summarize_run)
    monkeypatch.setattr("quacc.recipes.vasp.slabs.summarize_run", mock_summarize_run)
