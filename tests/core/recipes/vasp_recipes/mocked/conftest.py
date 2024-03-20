import warnings
from pathlib import Path

import pytest
from emmet.core.tasks import TaskDoc

FILE_DIR = Path(__file__).parent
PSEUDO_DIR = FILE_DIR / "fake_pseudos"

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    MOCK_TASKDOC = TaskDoc.from_directory(FILE_DIR / "mocked_vasp_run")


def mock_run(self, *args, **kwargs):
    from ase.io import write

    write(Path(self.directory, "CONTCAR"), self.atoms)


@pytest.fixture(autouse=True)
def patch_run(monkeypatch):
    from quacc.calculators.vasp.vasp import Vasp

    monkeypatch.setenv("VASP_PP_PATH", str(PSEUDO_DIR.as_posix()))
    monkeypatch.setattr(Vasp, "_run", mock_run)


def mock_read_results(self, *args, **kwargs):
    from ase.calculators.emt import EMT

    atoms = self.atoms
    atoms.calc = EMT()
    atoms.get_potential_energy()
    self.results.update(atoms.calc.results)


@pytest.fixture(autouse=True)
def patch_read_results(monkeypatch):
    from quacc.calculators.vasp.vasp import Vasp

    monkeypatch.setattr(Vasp, "read_results", mock_read_results)


def mock_taskdoc(dir_name, *args, **kwargs):
    from ase.io import read
    from monty.os.path import zpath

    from quacc.atoms.core import check_is_metal

    MOCK_TASKDOC.output.bandgap = (
        0.0 if check_is_metal(read(zpath(Path(dir_name, "CONTCAR")))) else 0.5
    )
    return MOCK_TASKDOC


@pytest.fixture(autouse=True)
def patch_taskdoc(monkeypatch):
    monkeypatch.setattr("quacc.schemas.vasp.TaskDoc.from_directory", mock_taskdoc)
