import pytest

dask = pytest.importorskip("dask")
from shutil import which

from dask.distributed import default_client

from quacc import SETTINGS

pytestmark = pytest.mark.skipif(
    which(str(SETTINGS.ESPRESSO_BIN_DIR / SETTINGS.ESPRESSO_BINARIES["pw"])) is None
    or which(str(SETTINGS.ESPRESSO_BIN_DIR / SETTINGS.ESPRESSO_BINARIES["ph"])) is None,
    reason="QE not installed",
)

from pathlib import Path

from ase.build import bulk

from quacc.recipes.espresso.phonons import grid_phonon_flow
from quacc.utils.files import copy_decompress_files

DATA_DIR = (
    Path(__file__).parent / ".." / "core" / "recipes" / "espresso_recipes" / "data"
)
client = default_client()


def test_phonon_grid_single(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "electrons": {"conv_thr": 1.0e-5},
        "control": {"pseudo_dir": tmp_path},
    }

    ph_loose = {
        "inputph": {"tr2_ph": 1e-6, "qplot": False, "ldisp": False, "lqdir": True}
    }

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {
        "relax_job": {
            "input_data": input_data,
            "pseudopotentials": pseudopotentials,
            "kspacing": 0.5,
        },
        # Doing a mistake in qpts here for testing purposes
        "ph_job": {"input_data": ph_loose, "qpts": [(0.1, 0, 0, 1)]},
    }

    future = grid_phonon_flow(atoms, job_params=job_params)
    grid_results = client.compute(future).result()
    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in grid_results["results"][1]


def test_phonon_grid_single_gamma(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "electrons": {"conv_thr": 1.0e-5},
        "control": {"pseudo_dir": tmp_path},
    }

    ph_loose = {
        "inputph": {"tr2_ph": 1e-6, "qplot": False, "ldisp": False, "lqdir": False}
    }

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {
        "relax_job": {
            "input_data": input_data,
            "pseudopotentials": pseudopotentials,
            "kspacing": 0.5,
        },
        "ph_job": {"input_data": ph_loose, "qpts": (0.0, 0, 0)},
    }

    future = grid_phonon_flow(atoms, job_params=job_params)
    grid_results = client.compute(future).result()
    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in grid_results["results"][1]


def test_phonon_grid_qplot(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "electrons": {"conv_thr": 1.0e-5},
        "control": {"pseudo_dir": tmp_path},
    }

    ph_loose = {
        "inputph": {"tr2_ph": 1e-6, "qplot": True, "ldisp": True, "lqdir": True}
    }

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {
        "relax_job": {
            "input_data": input_data,
            "pseudopotentials": pseudopotentials,
            "kspacing": 0.5,
        },
        "ph_job": {"input_data": ph_loose, "qpts": [(0.1, 0, 0, 1), (0.2, 0, 0, 1)]},
    }

    future = grid_phonon_flow(atoms, job_params=job_params)
    grid_results = client.compute(future).result()

    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in grid_results["results"][1]
        assert key in grid_results["results"][2]


def test_phonon_grid_disp(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "electrons": {"conv_thr": 1.0e-1},
        "control": {"pseudo_dir": tmp_path},
    }

    ph_loose = {
        "inputph": {
            "tr2_ph": 1e-2,
            "qplot": False,
            "ldisp": True,
            "lqdir": True,
            "nq1": 2,
            "nq2": 2,
            "nq3": 2,
        }
    }

    pseudopotentials = {"Si": "Si.upf"}

    job_params = {
        "relax_job": {
            "input_data": input_data,
            "pseudopotentials": pseudopotentials,
            "kspacing": 0.5,
        },
        "ph_job": {"input_data": ph_loose},
    }

    future = grid_phonon_flow(atoms, job_params=job_params)
    grid_results = client.compute(future).result()

    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in grid_results["results"][1]


def test_phonon_grid_v2(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files([DATA_DIR / "Li.upf.gz"], tmp_path)

    atoms = bulk("Li", "bcc", orthorhombic=True)

    input_data = {
        "electrons": {"conv_thr": 1.0e-5},
        "control": {"pseudo_dir": tmp_path},
    }

    ph_loose = {"inputph": {"tr2_ph": 1e-6, "lqdir": True}}

    pseudopotentials = {"Li": "Li.upf"}

    job_params = {
        "relax_job": {
            "input_data": input_data,
            "pseudopotentials": pseudopotentials,
            "kspacing": 0.5,
        },
        "ph_job": {"input_data": ph_loose, "qpts": (0.0, 0.0, 0.0)},
    }

    future = grid_phonon_flow(atoms, job_params=job_params, nblocks=3)
    grid_results = client.compute(future).result()

    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in grid_results["results"][1]
