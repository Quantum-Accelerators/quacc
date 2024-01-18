
import pytest
parsl = pytest.importorskip("parsl")
from shutil import which

pytestmark = pytest.mark.skipif(
    which("pw.x") is None or which("ph.x") is None, reason="QE not installed"
)
import contextlib
from pathlib import Path

from ase.build import bulk

from quacc import SETTINGS
from quacc.recipes.espresso.phonons import grid_phonon_flow
from quacc.utils.files import copy_decompress_files

DEFAULT_SETTINGS = SETTINGS.model_copy()
DATA_DIR = Path(__file__).parent /".."/"core"/"recipes"/"espresso_recipes"/ "data"



def setup_module():
    with contextlib.suppress(Exception):
        parsl.load()


def teardown_module():
    parsl.clear()
def test_phonon_grid(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "electrons": {"conv_thr": 1.0e-5},
        "control": {"pseudo_dir": tmp_path},
    }

    ph_loose = {"inputph": {"tr2_ph": 1e-6, "qplot": True, "ldisp": True}}

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {
        "relax_job": {
            "input_data": input_data,
            "pseudopotentials": pseudopotentials,
            "kspacing": 0.5,
        },
        "ph_job": {"input_data": ph_loose, "qpts": [(0, 0, 0, 1), (0.5, 0.0, 0.0, 1)]},
    }

    grid_results = grid_phonon_flow(atoms, job_params=job_params)

    assert grid_results.result()