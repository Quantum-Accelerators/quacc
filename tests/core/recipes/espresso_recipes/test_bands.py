from pathlib import Path
from shutil import which

import pytest
from ase.build import bulk
from numpy.testing import assert_allclose

from quacc.recipes.espresso.bands import bands_job
from quacc.utils.files import copy_decompress_files, copy_decompress_tree

pytestmark = pytest.mark.skipif(
    which("pw.x") is None or which("bands.x") is None, reason="QE not installed"
)

DATA_DIR = Path(__file__).parent / "data"


def test_bands_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_tree({DATA_DIR / "dos_test/": "pwscf.save/*.gz"}, tmp_path)
    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)
    atoms = bulk("Si")
    pseudopotentials = {"Si": "Si.upf"}
    job_params = {
        "bands": {
            "input_data": {"control": {"pseudo_dir": tmp_path}},
            "pseudopotentials": pseudopotentials,
        },
        "bands_pp": {},
        "fermi_surface": {},
    }

    output = bands_job(atoms, prev_dir=tmp_path, job_params=job_params)
    assert (
        output["bands"]["parameters"]["input_data"]["control"]["calculation"] == "bands"
    )

    assert_allclose(
        output["bands"]["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
    )

    assert output["bands"]["results"]["nbands"] == 4

    assert_allclose(
        output["bands_pp"]["atoms"].get_positions(), atoms.get_positions(), atol=1.0e-4
    )

    assert output["bands_pp"]["name"] == "bands.x post-processing"
