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
from numpy.testing import assert_allclose, assert_array_equal

from quacc.recipes.espresso.core import static_job
from quacc.recipes.espresso.phonons import phonon_job
from quacc.utils.files import copy_decompress_files

DEFAULT_SETTINGS = SETTINGS.model_copy()
DATA_DIR = Path(__file__).parent / "data"


def test_phonon_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Li")

    copy_decompress_files([DATA_DIR / "Li.upf.gz"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    input_data = {
        "control": {"calculation": "scf"},
        "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.02},
        "electrons": {"mixing_mode": "TF", "mixing_beta": 0.7, "conv_thr": 1.0e-6},
    }

    ph_loose = {"inputph": {"tr2_ph": 1e-8}}

    pseudopotentials = {"Li": "Li.upf"}

    pw_results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    ph_results = phonon_job(pw_results["dir_name"], input_data=ph_loose)

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


def test_phonon_job_list_to_do(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Li")

    copy_decompress_files([DATA_DIR / "Li.upf.gz"], tmp_path)

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
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    qpts = [(0, 0, 0, 1), (1 / 3, 0, 0, 1), (1 / 2, 0, 0, 1)]

    nat_todo = [1]

    ph_results = phonon_job(
        pw_results["dir_name"],
        input_data=ph_loose,
        qpts=qpts,
        nat_todo_indices=nat_todo,
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
