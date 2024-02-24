from pathlib import Path

import pytest
from ase.calculators.emt import EMT
from ase.io import write
from emmet.core.tasks import TaskDoc

FILE_DIR = Path(__file__).parent
PSEUDO_DIR = FILE_DIR / "fake_pseudos"

MOCK_TASKDOC = TaskDoc.from_directory(FILE_DIR / "mocked_vasp_run")


def mock_run(self, *args, **kwargs):
    write(Path(self.directory) / "CONTCAR", self.atoms)


@pytest.fixture(autouse=True)
def patch_run(monkeypatch):
    from quacc.calculators.vasp.vasp import Vasp

    monkeypatch.setenv("VASP_PP_PATH", str(PSEUDO_DIR))
    monkeypatch.setattr(Vasp, "_run", mock_run)


def mock_read_results(self, *args, **kwargs):
    atoms = self.atoms
    atoms.calc = EMT()
    atoms.get_potential_energy()
    self.results.update(atoms.calc.results)


@pytest.fixture(autouse=True)
def patch_read_results(monkeypatch):
    from quacc.calculators.vasp.vasp import Vasp

    monkeypatch.setattr(Vasp, "read_results", mock_read_results)


def mock_taskdoc(*args, **kwargs):
    return MOCK_TASKDOC


@pytest.fixture(autouse=True)
def patch_taskdoc(monkeypatch):
    monkeypatch.setattr("quacc.schemas.vasp.TaskDoc.from_directory", mock_taskdoc)
