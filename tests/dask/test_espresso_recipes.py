from __future__ import annotations

import pytest

dask = pytest.importorskip("dask")
pytest.importorskip("distributed")

import gzip
from shutil import which

from dask.distributed import default_client

from quacc import get_settings

settings = get_settings()
pytestmark = pytest.mark.skipif(
    which(str(settings.ESPRESSO_BIN_DIR / settings.ESPRESSO_BINARIES["pw"])) is None
    or which(str(settings.ESPRESSO_BIN_DIR / settings.ESPRESSO_BINARIES["ph"])) is None,
    reason="QE not installed",
)

from pathlib import Path

from ase.build import bulk

from quacc import subflow
from quacc.recipes.espresso.core import post_processing_job, relax_job, static_job
from quacc.recipes.espresso.phonons import grid_phonon_flow
from quacc.utils.files import copy_decompress_files

DATA_DIR = (
    Path(__file__).parent / ".." / "core" / "recipes" / "espresso_recipes" / "data"
)
client = default_client()


def test_phonon_grid_single(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR, ["Si.upf.gz"], tmp_path)

    ph_loose = {
        "inputph": {"tr2_ph": 1e-6, "qplot": False, "ldisp": False, "lqdir": True}
    }

    relax_output = client.compute(
        relax_job(
            bulk("Si"),
            input_data={"control": {"pseudo_dir": tmp_path}},
            pseudopotentials={"Si": "Si.upf"},
            kspacing=0.5,
        )
    ).result()
    # Doing a mistake in qpts here for testing purposes
    job_params = {"ph_job": {"input_data": ph_loose, "qpts": [(0.1, 0, 0, 1)]}}

    future = grid_phonon_flow(
        prev_outdir=relax_output["dir_name"], job_params=job_params
    )
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

    ph_loose = {
        "inputph": {"tr2_ph": 1e-6, "qplot": False, "ldisp": False, "lqdir": False}
    }

    relax_output = client.compute(
        relax_job(
            bulk("Si"),
            input_data={"control": {"pseudo_dir": tmp_path}},
            pseudopotentials={"Si": "Si.upf"},
            kspacing=0.5,
        )
    ).result()
    job_params = {"ph_job": {"input_data": ph_loose, "qpts": (0.0, 0, 0)}}

    future = grid_phonon_flow(
        prev_outdir=relax_output["dir_name"], job_params=job_params
    )
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

    ph_loose = {
        "inputph": {"tr2_ph": 1e-6, "qplot": True, "ldisp": True, "lqdir": True}
    }

    relax_output = client.compute(
        relax_job(
            bulk("Si"),
            input_data={"control": {"pseudo_dir": tmp_path}},
            pseudopotentials={"Si": "Si.upf"},
            kspacing=0.5,
        )
    ).result()
    job_params = {
        "ph_job": {"input_data": ph_loose, "qpts": [(0.1, 0, 0, 1), (0.2, 0, 0, 1)]}
    }

    future = grid_phonon_flow(
        prev_outdir=relax_output["dir_name"], job_params=job_params
    )
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

    relax_output = client.compute(
        relax_job(
            bulk("Si"),
            input_data={
                "electrons": {"conv_thr": 1.0e-1},
                "control": {"pseudo_dir": tmp_path},
            },
            pseudopotentials={"Si": "Si.upf"},
            kspacing=0.5,
        )
    ).result()
    job_params = {"ph_job": {"input_data": ph_loose}}

    future = grid_phonon_flow(
        prev_outdir=relax_output["dir_name"], job_params=job_params
    )
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

    ph_loose = {"inputph": {"tr2_ph": 1e-6, "lqdir": True}}

    relax_output = client.compute(
        relax_job(
            bulk("Li", "bcc", orthorhombic=True),
            input_data={"control": {"pseudo_dir": tmp_path}},
            pseudopotentials={"Li": "Li.upf"},
            kspacing=0.5,
        )
    ).result()
    job_params = {"ph_job": {"input_data": ph_loose, "qpts": (0.0, 0.0, 0.0)}}

    future = grid_phonon_flow(
        prev_outdir=relax_output["dir_name"], job_params=job_params, nblocks=3
    )
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
        return [
            post_processing_job(prev_outdir=results, input_data={"plot_num": plot_num})
            for plot_num in [0, 1, 2, 4, 8, 123, 3]
        ]

    future = pp_subflow(static_results["dir_name"])
    pp_results = client.compute(future).result()

    for pp_result in pp_results:
        assert Path(pp_result["dir_name"], "pseudo_charge_density.cube.gz").is_file()

        assert (
            pp_result["parameters"]["input_data"]["inputpp"]["outdir"]
            == static_results["dir_name"]
        )
