from __future__ import annotations

from pathlib import Path

import pytest
from ase.build import bulk
from numpy.testing import assert_allclose

from quacc import SETTINGS
from quacc.recipes.espresso.dos import dos_flow
from quacc.utils.files import copy_decompress_files

DEFAULT_SETTINGS = SETTINGS.model_copy()
DATA_DIR = Path(__file__).parent / "data"


def test_dos_flow(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)
    SETTINGS.ESPRESSO_PSEUDO = tmp_path
    atoms = bulk("Si")
    input_data = {
        "control": {"calculation": "scf"},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {
        "static_job": {"input_data": input_data, "pseudopotentials": pseudopotentials},
        "non_scf_job": {"kspacing": 0.05},
    }

    output = dos_flow(atoms, job_params=job_params)
    assert_allclose(
        output["static_job"]["atoms"].get_positions(),
        atoms.get_positions(),
        atol=1.0e-4,
    )
    assert (
        output["static_job"]["parameters"]["input_data"]["control"]["calculation"]
        == "scf"
    )
    assert (
        output["static_job"]["parameters"]["input_data"]["electrons"]["mixing_mode"]
        == "TF"
    )
    assert_allclose(
        output["non_scf_job"]["atoms"].get_positions(),
        atoms.get_positions(),
        atol=1.0e-4,
    )
    assert (
        output["non_scf_job"]["parameters"]["input_data"]["control"]["calculation"]
        == "nscf"
    )
    assert (
        output["non_scf_job"]["parameters"]["input_data"]["electrons"]["mixing_mode"]
        == "TF"
    )
