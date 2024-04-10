from __future__ import annotations

from pathlib import Path

import pytest

FILE_DIR = Path(__file__).parent
ORCA_DIR = Path(FILE_DIR, "orca_run")


def mock_execute(self, directory, *args, **kwargs):
    from shutil import copy

    copy(ORCA_DIR / "orca.out.gz", Path(directory, "orca.out.gz"))


@pytest.fixture(autouse=True)
def patch_execute(monkeypatch):
    from ase.calculators.orca import OrcaTemplate

    monkeypatch.setattr(OrcaTemplate, "execute", mock_execute)


def mock_read_results(self, directory, *args, **kwargs):
    from ase.calculators.lj import LennardJones
    from ase.io import write
    from ase.io.orca import read_geom_orcainp

    atoms = read_geom_orcainp(directory / "orca.inp")
    write(directory / "orca.xyz", atoms)
    atoms.calc = LennardJones()
    atoms.get_potential_energy()
    return atoms.calc.results


@pytest.fixture(autouse=True)
def patch_read_results(monkeypatch):
    from ase.calculators.orca import OrcaTemplate

    monkeypatch.setattr(OrcaTemplate, "read_results", mock_read_results)
