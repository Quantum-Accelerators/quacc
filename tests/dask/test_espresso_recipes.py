from __future__ import annotations

import pytest

dask = pytest.importorskip("dask")
pytest.importorskip("distributed")

import gzip
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

from quacc import subflow
from quacc.recipes.espresso.bands import bands_pw_job
from quacc.recipes.espresso.core import post_processing_job, static_job
from quacc.recipes.espresso.phonons import grid_phonon_flow, phonon_job, postahc_job
from quacc.utils.files import copy_decompress_files

DATA_DIR = (
    Path(__file__).parent / ".." / "core" / "recipes" / "espresso_recipes" / "data"
)
client = default_client()


def test_phonon_grid_single(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR, ["Si.upf.gz"], tmp_path)

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

    outfile = Path(grid_results["dir_name"], "ph.out.gz")

    assert outfile.exists()

    with gzip.open(outfile, "r") as fd:
        lines = str(fd.read())

    assert "Self-consistent Calculation" not in lines


def test_phonon_grid_single_gamma(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR, ["Si.upf.gz"], tmp_path)

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
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR, ["Si.upf.gz"], tmp_path)

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
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR, ["Si.upf.gz"], tmp_path)

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
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

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


def test_pp_concurrent_inplace(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR, ["Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "system": {"occupations": "smearing", "smearing": "gaussian", "degauss": 0.005},
        "electrons": {"mixing_mode": "plain", "mixing_beta": 0.6, "conv_thr": 1.0e-6},
        "control": {"pseudo_dir": tmp_path},
    }

    pseudopotentials = {"Si": "Si.upf"}

    static_results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    static_results = client.compute(static_results).result()

    @subflow
    def pp_subflow(results):
        pp_results = []

        for plot_num in [0, 1, 2, 4, 8, 123, 3]:
            pp_results.append(
                post_processing_job(
                    prev_outdir=results, input_data={"plot_num": plot_num}
                )
            )

        return pp_results

    future = pp_subflow(static_results["dir_name"])
    pp_results = client.compute(future).result()

    for pp_result in pp_results:
        assert Path(pp_result["dir_name"], "pseudo_charge_density.cube.gz").is_file()

        assert (
            pp_result["parameters"]["input_data"]["inputpp"]["outdir"]
            == static_results["dir_name"]
        )


def test_bands_pw_concurrent_inplace(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR, ["Si.upf.gz"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "system": {"occupations": "smearing", "smearing": "gaussian", "degauss": 0.005},
        "electrons": {"mixing_mode": "plain", "mixing_beta": 0.6, "conv_thr": 1.0e-6},
        "control": {"pseudo_dir": tmp_path},
    }

    pseudopotentials = {"Si": "Si.upf"}

    static_results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    @subflow
    def bands_subflow(copy_files):
        bands_results = []

        for line_density in [10, 20, 30]:
            bands_results.append(
                bands_pw_job(atoms, copy_files, line_density=line_density)
            )

        return bands_results

    future = bands_subflow(static_results["dir_name"])
    bands_results = client.compute(future).result()

    for bands_result in bands_results:
        assert bands_result["name"] == "pw.x bands"
        assert bands_result["results"]["nbands"] == 4


def test_postahc_concurrent(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR, ["Si.upf.gz"], tmp_path)
    atoms = bulk("Si")

    input_data = {
        "control": {"calculation": "scf", "pseudo_dir": tmp_path},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }
    pseudopotentials = {"Si": "Si.upf"}

    static_results = static_job(
        atoms,
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kspacing=0.5,
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    ph_input_data = {
        "inputph": {"tr2_ph": 1e-12, "electron_phonon": "ahc", "ahc_nbnd": 4}
    }
    phonon_results = phonon_job(
        static_results["dir_name"],
        input_data=ph_input_data,
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    @subflow
    def postahc_subflow(copy_files):
        postahc_results = []

        for eta in [0.01, 0.02]:
            postahc_results.append(
                postahc_job(
                    copy_files,
                    input_data={"eta": eta},
                    parallel_info=ESPRESSO_PARALLEL_INFO,
                )
            )

        return postahc_results

    future = postahc_subflow(phonon_results["dir_name"])
    postahc_results = client.compute(future).result()

    for postahc_result in postahc_results:
        assert Path(postahc_result["dir_name"], "postahc.out.gz").exists()
        assert Path(postahc_result["dir_name"], "selfen_real.dat.gz").exists()
