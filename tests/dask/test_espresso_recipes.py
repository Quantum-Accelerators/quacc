import pytest
from shutil import which

dask = pytest.importorskip("dask")

pytestmark = pytest.mark.skipif(
    which("pw.x") is None or which("dos.x") is None, reason="QE not installed"
)

from dask.distributed import default_client

from pathlib import Path

import pytest
from ase.build import bulk

from quacc import SETTINGS
from quacc.recipes.espresso.dos import dos_flow
from quacc.utils.files import copy_decompress_files

client = default_client()

DEFAULT_SETTINGS = SETTINGS.model_copy()
DATA_DIR = (
    Path(__file__).parent / ".." / "core" / "recipes" / "espresso_recipes" / "data"
)


def test_dos_flow(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)
    atoms = bulk("Si")
    input_data = {
        "control": {"calculation": "scf", "pseudo_dir": tmp_path},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {
        "static_job": {"input_data": input_data, "pseudopotentials": pseudopotentials},
        "non_scf_job": {"kspacing": 0.05},
    }

    future = dos_flow(atoms, job_params=job_params)
    assert client.compute(future).result()
