from pathlib import Path
from shutil import copy

import pytest
from ase.calculators.lj import LennardJones
from ase.calculators.orca import OrcaTemplate
from ase.io import write
from ase.io.orca import read_geom_orcainp

FILE_DIR = Path(__file__).parent
ORCA_DIR = Path(FILE_DIR, "orca_run")


def mock_execute(self, *args, **kwargs):
    copy(ORCA_DIR / "orca.out.gz", "orca.out.gz")


@pytest.fixture(autouse=True)
def patch_execute(monkeypatch):
    monkeypatch.setattr(OrcaTemplate, "execute", mock_execute)


def mock_read_results(self, directory, *args, **kwargs):
    atoms = read_geom_orcainp(directory / "orca.inp")
    write(directory / "orca.xyz", atoms)
    atoms.calc = LennardJones()
    atoms.get_potential_energy()
    return atoms.calc.results


@pytest.fixture(autouse=True)
def patch_read_results(monkeypatch):
    monkeypatch.setattr(OrcaTemplate, "read_results", mock_read_results)
