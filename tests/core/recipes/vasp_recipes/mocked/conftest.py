from pathlib import Path

import pytest
from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.vasp.vasp import Vasp
from ase.io import write

FILE_DIR = Path(__file__).parent
PSEUDO_DIR = FILE_DIR / "fake_pseudos"


def mock_execute(self, *args, **kwargs):
    write(Path(self.directory) / "CONTCAR", self.atoms)


@pytest.fixture(autouse=True)
def patch_execute(monkeypatch):
    monkeypatch.setenv("VASP_PP_PATH", str(PSEUDO_DIR))
    monkeypatch.setattr(Vasp, "_run", mock_execute)


def mock_read_results(self, *args, **kwargs):
    atoms = self.atoms
    atoms.calc = EMT()
    atoms.get_potential_energy()
    self.results.update(atoms.calc.results)


@pytest.fixture(autouse=True)
def patch_read_results(monkeypatch):
    monkeypatch.setattr(Vasp, "read_results", mock_read_results)


def mock_summarize_run(atoms, **kwargs):
    from quacc.atoms.core import check_is_metal
    from quacc.schemas.ase import summarize_run as calc_summarize_run

    # Instead of running the VASP-specific summarize_run(), we mock it with the
    # general calculator schema which does not require VASP files to be
    # in the working directory and will work with pytest.

    move_magmoms = kwargs.get("move_magmoms", True)
    additional_fields = kwargs.get("additional_fields")
    output = calc_summarize_run(
        atoms, Atoms(), move_magmoms=move_magmoms, additional_fields=additional_fields
    )
    output["output"] = {
        "energy": -1.0,
        "bandgap": 0.0 if check_is_metal(atoms) else 0.5,
    }
    return output


@pytest.fixture(autouse=True)
def patch_summarize_run(monkeypatch):
    # Monkeypatch the summarize_run() function so that we aren't relying on real
    # VASP files to be in the working directory during the test.
    monkeypatch.setattr(
        "quacc.recipes.vasp._base.vasp_summarize_run", mock_summarize_run
    )
