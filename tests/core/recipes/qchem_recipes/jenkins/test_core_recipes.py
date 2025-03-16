from __future__ import annotations

from shutil import which

import pytest

from quacc import get_settings

pytestmark = pytest.mark.skipif(
    which(str(get_settings().QCHEM_CMD)) is None, reason="QChem not installed"
)
from pathlib import Path

from ase.build import molecule

from quacc.recipes.qchem.core import relax_job, static_job

FILE_DIR = Path(__file__).parent
QCHEM_DIR = FILE_DIR / "qchem_examples"


def test_static_job_simple(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    output = static_job(
        molecule("O2"), spin_multiplicity=3, method="pbe", basis="def2-svp"
    )
    assert output["molecule_metadata"]["formula_alphabetical"] == "O2"
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 3
    assert output["results"]["energy"] < 0
    assert output["atoms"] == molecule("O2")


def test_relax_job_simple(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    output = relax_job(
        molecule("O2"), spin_multiplicity=3, method="pbe", basis="def2-svp"
    )

    assert output["molecule_metadata"]["formula_alphabetical"] == "O2"
    assert output["parameters"]["charge"] == 0
    assert output["parameters"]["spin_multiplicity"] == 3
    assert output["results"]["energy"] < 0
    assert len(output["results"]["taskdoc"]) > 1
    assert output["atoms"] != molecule("O2")
