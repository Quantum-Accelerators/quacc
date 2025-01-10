from __future__ import annotations

import warnings
from pathlib import Path

import pytest
from emmet.core.tasks import TaskDoc

FILE_DIR = Path(__file__).parent
PSEUDO_DIR = FILE_DIR / "fake_pseudos"

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    MOCK_METALLIC_TASKDOC = TaskDoc.from_directory(
        FILE_DIR / "mocked_vasp_runs" / "metallic"
    )
    MOCK_NONMETALLIC_TASKDOC = TaskDoc.from_directory(
        FILE_DIR / "mocked_vasp_runs" / "nonmetallic"
    )


def mock_run(self, *args, **kwargs):
    from ase.io import write

    write(Path(self.directory, "CONTCAR"), self.atoms)
    return None, None


@pytest.fixture(autouse=True)
def patch_run(monkeypatch):
    from quacc.calculators.vasp.vasp import Vasp

    monkeypatch.setenv("VASP_PP_PATH", str(PSEUDO_DIR))
    monkeypatch.setattr(Vasp, "_run", mock_run)


def mock_read_results(self, *args, **kwargs):
    from ase.calculators.lj import LennardJones

    atoms = self.atoms
    atoms.calc = LennardJones()
    atoms.get_potential_energy()
    self.results.update(atoms.calc.results)


@pytest.fixture(autouse=True)
def patch_read_results(monkeypatch):
    from quacc.calculators.vasp.vasp import Vasp

    monkeypatch.setattr(Vasp, "read_results", mock_read_results)


def mock_metallic_taskdoc(dir_name, *args, **kwargs):
    return MOCK_METALLIC_TASKDOC


@pytest.fixture
def patch_metallic_taskdoc(monkeypatch):
    monkeypatch.setattr(
        "quacc.schemas.vasp.TaskDoc.from_directory", mock_metallic_taskdoc
    )

    class DummyBandStructure:
        def __init__(self, *args, **kwargs):
            pass

        @staticmethod
        def is_metal():
            return True

        @staticmethod
        def get_band_gap():
            return {"energy": MOCK_METALLIC_TASKDOC.output.bandgap}

    monkeypatch.setattr(
        "pymatgen.io.vasp.Vasprun.get_band_structure", DummyBandStructure
    )


def mock_nonmetallic_taskdoc(dir_name, *args, **kwargs):
    return MOCK_NONMETALLIC_TASKDOC


@pytest.fixture
def patch_nonmetallic_taskdoc(monkeypatch):
    monkeypatch.setattr(
        "quacc.schemas.vasp.TaskDoc.from_directory", mock_nonmetallic_taskdoc
    )

    class DummyBandStructure:
        def __init__(self, *args, **kwargs):
            pass

        @staticmethod
        def is_metal():
            return False

        @staticmethod
        def get_band_gap():
            return {"energy": MOCK_NONMETALLIC_TASKDOC.output.bandgap}

    monkeypatch.setattr(
        "pymatgen.io.vasp.Vasprun.get_band_structure", DummyBandStructure
    )
