from shutil import which

import pytest

from quacc import SETTINGS

pytestmark = pytest.mark.skipif(
    which(str(SETTINGS.ESPRESSO_BIN_DIR / SETTINGS.ESPRESSO_BINARIES["pw"])) is None
    or which(str(SETTINGS.ESPRESSO_BIN_DIR / SETTINGS.ESPRESSO_BINARIES["ph"])) is None,
    reason="QE not installed",
)

from pathlib import Path

from ase.build import bulk
from ase.io.espresso import read_fortran_namelist
from monty.shutil import decompress_file
from numpy.testing import assert_allclose, assert_array_equal

from quacc.recipes.espresso.core import static_job
from quacc.recipes.espresso.phonons import (
    matdyn_job,
    phonon_dos_flow,
    phonon_job,
    q2r_job,
)
from quacc.utils.files import copy_decompress_files

DEFAULT_SETTINGS = SETTINGS.model_copy()
DATA_DIR = Path(__file__).parent / "data"


def test_phonon_job(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    atoms = bulk("Li")

    copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    input_data = {
        "control": {"calculation": "scf"},
        "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    ph_loose = {"inputph": {"tr2_ph": 1e-8}}

    pseudopotentials = {"Li": "Li.upf"}

    pw_results = static_job(
        atoms,
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kspacing=0.5,
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    ph_results = phonon_job(
        pw_results["dir_name"], input_data=ph_loose, parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    assert_allclose(
        ph_results["results"][1]["atoms"].get_positions(),
        atoms.get_positions(),
        atol=1.0e-4,
    )
    # ph.x cell param are not defined to a very high level of accuracy,
    # atol = 1.0e-3 is needed here...
    assert_allclose(
        ph_results["results"][1]["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3
    )
    assert_array_equal(
        ph_results["results"][1]["atoms"].get_chemical_symbols(),
        atoms.get_chemical_symbols(),
    )

    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in ph_results["results"][1]

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO


def test_phonon_job_list_to_do(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    atoms = bulk("Li")

    copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    input_data = {
        "control": {"calculation": "scf"},
        "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }
    ph_loose = {
        "inputph": {"tr2_ph": 1e-8, "qplot": True, "nat_todo": 1, "ldisp": True}
    }
    pseudopotentials = {"Li": "Li.upf"}

    pw_results = static_job(
        atoms,
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kspacing=0.5,
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    qpts = [(0, 0, 0, 1), (1 / 3, 0, 0, 1), (1 / 2, 0, 0, 1)]

    nat_todo = [1]

    ph_results = phonon_job(
        pw_results["dir_name"],
        input_data=ph_loose,
        qpts=qpts,
        nat_todo_indices=nat_todo,
        parallel_info=ESPRESSO_PARALLEL_INFO,
    )

    assert_allclose(
        ph_results["results"][1]["atoms"].get_positions(),
        atoms.get_positions(),
        atol=1.0e-4,
    )
    # ph.x cell param are not defined to a very high level of accuracy,
    # atol = 1.0e-3 is needed here...
    assert_allclose(
        ph_results["results"][1]["atoms"].get_cell(), atoms.get_cell(), atol=1.0e-3
    )
    assert_array_equal(
        ph_results["results"][1]["atoms"].get_chemical_symbols(),
        atoms.get_chemical_symbols(),
    )

    sections = [
        "atoms",
        "eqpoints",
        "freqs",
        "kpoints",
        "mode_symmetries",
        "representations",
    ]

    for key in sections:
        assert key in ph_results["results"][1]

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO


def test_q2r_job(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR / "q2r_test", "matdyn", tmp_path)

    additional_cards = ["1 1 1", "1", "matdyn"]

    q2r_results = q2r_job(
        tmp_path, additional_cards=additional_cards, parallel_info=ESPRESSO_PARALLEL_INFO
    )

    assert Path(q2r_results["dir_name"], "q2r.fc.gz").exists()

    decompress_file(Path(q2r_results["dir_name"], "q2r.in.gz"))

    with Path(q2r_results["dir_name"], "q2r.in").open() as f:
        recycled_input = read_fortran_namelist(f)

    assert Path(recycled_input[0]["input"].pop("flfrc")).is_absolute()

    assert recycled_input[0]["input"] == {"fildyn": "matdyn"}
    assert recycled_input[1] == additional_cards + ["EOF"]


def test_matdyn_job(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    copy_decompress_files(DATA_DIR / "matdyn_test", "q2r.fc", tmp_path)

    input_data = {"input": {"dos": True, "nk1": 4, "nk2": 4, "nk3": 4}}
    matdyn_results = matdyn_job(
        tmp_path, input_data=input_data, parallel_info=ESPRESSO_PARALLEL_INFO
    )

    assert Path(matdyn_results["dir_name"], "q2r.fc.gz").exists()
    assert Path(matdyn_results["dir_name"], "matdyn.dos.gz").exists()
    assert Path(matdyn_results["dir_name"], "matdyn.freq.gz").exists()
    assert Path(matdyn_results["dir_name"], "matdyn.modes.gz").exists()
    assert matdyn_results["results"]["matdyn_dos"]["phonon_dos"].shape == (561, 3)

    decompress_file(Path(matdyn_results["dir_name"], "matdyn.in.gz"))

    with Path(matdyn_results["dir_name"], "matdyn.in").open() as f:
        recycled_input = read_fortran_namelist(f)

    assert Path(recycled_input[0]["input"].pop("flfrc")).is_absolute()

    assert recycled_input[0]["input"] == input_data["input"]


def test_phonon_dos_flow(tmp_path, monkeypatch, ESPRESSO_PARALLEL_INFO):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    atoms = bulk("Li")

    copy_decompress_files(DATA_DIR, ["Li.upf.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    input_data = {
        "control": {"calculation": "scf"},
        "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    pseudopotentials = {"Li": "Li.upf"}

    job_params = {
        "relax_job": {
            "pseudopotentials": pseudopotentials,
            "input_data": input_data,
            "kspacing": 1.0,
            "parallel_info": ESPRESSO_PARALLEL_INFO,
        },
        "phonon_job": {"parallel_info": ESPRESSO_PARALLEL_INFO},
        "q2r_job": {"parallel_info": ESPRESSO_PARALLEL_INFO},
        "matdyn_job": {"parallel_info": ESPRESSO_PARALLEL_INFO},
    }

    phonon_dos_results = phonon_dos_flow(atoms, job_params=job_params)

    SETTINGS.ESPRESSO_PSEUDO = DEFAULT_SETTINGS.ESPRESSO_PSEUDO
