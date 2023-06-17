import os
from pathlib import Path
from shutil import copy, rmtree

import numpy as np
import pytest
from ase.calculators.calculator import FileIOCalculator

FILE_DIR = Path(__file__).resolve().parent
QCHEM_DIR = os.path.join(FILE_DIR, "qchem_examples")
TMP_DIR = os.path.join(FILE_DIR, "tmp_dir")


def mock_execute(self, **kwargs):
    if not os.path.exists(os.path.join(FILE_DIR, "mol.qout.gz")):
        copy(
            os.path.join(QCHEM_DIR, "mol.qout.basic"), os.path.join(TMP_DIR, "mol.qout")
        )
        copy(os.path.join(QCHEM_DIR, "131.0.basic"), os.path.join(TMP_DIR, "131.0"))
    else:
        copy(
            os.path.join(QCHEM_DIR, "mol.qout.intermediate"),
            os.path.join(TMP_DIR, "mol.qout"),
        )
        copy(
            os.path.join(QCHEM_DIR, "131.0.intermediate"),
            os.path.join(TMP_DIR, "131.0"),
        )


@pytest.fixture(autouse=True)
def patch_execute(monkeypatch):
    # Monkeypatch the .execute() method of the FileIOCalculator object so
    # we aren't running the actual calculation during testing.
    monkeypatch.setattr(FileIOCalculator, "execute", mock_execute)
