from __future__ import annotations

from pathlib import Path
from shutil import which

import numpy as np
import pytest
from ase.build import bulk
from ase.io.espresso import construct_namelist

from quacc import SETTINGS
from quacc.recipes.espresso.core import ph_job, static_job
from quacc.utils.files import copy_decompress_files

pytestmark = pytest.mark.skipif(which("pw.x") is None, reason="QE not installed")

DEFAULT_SETTINGS = SETTINGS.model_copy()


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Si.upf"], tmp_path)

    atoms = bulk("Si")

    input_data = {
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.005,
        "mixing_mode": "plain",
        "mixing_beta": 0.6,
        "pseudo_dir": tmp_path,
        "conv_thr": 1.0e-6,
    }

    pseudopotentials = {"Si": "Si.upf"}

    results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.5
    )

    input_data = dict(construct_namelist(input_data))

    assert np.allclose(results["atoms"].positions, atoms.positions, atol=1.0e-4)

    assert np.allclose(results["atoms"].cell, atoms.cell, atol=1.0e-3)
    assert (results["atoms"].symbols == atoms.symbols).all()

    new_input_data = results["parameters"]["input_data"]

    assert new_input_data["system"]["degauss"] == 0.005
    assert new_input_data["system"]["occupations"] == "smearing"
    assert new_input_data["electrons"]["conv_thr"] == 1.0e-6
    assert new_input_data["control"]["calculation"] == "scf"


def test_ph_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Li")

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Li.upf"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = tmp_path

    input_data = {
        "calculation": "scf",
        "occupations": "smearing",
        "smearing": "cold",
        "degauss": 0.02,
        "mixing_mode": "TF",
        "mixing_beta": 0.7,
        "conv_thr": 1.0e-6,
    }

    ph_loose = {"tr2_ph": 1e-8}

    pseudopotentials = {"Li": "Li.upf"}

    pw_results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.25
    )

    ph_results = ph_job(input_data=ph_loose, copy_files=pw_results["dir_name"])

    assert (0, 0, 0) in ph_results["results"]
    assert np.allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].positions,
        atoms.positions,
        atol=1.0e-4,
    )
    # ph.x cell param are not defined to a very high level of accuracy,
    # atol = 1.0e-3 is needed here...
    assert np.allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].cell, atoms.cell, atol=1.0e-3
    )
    assert (ph_results["results"][(0, 0, 0)]["atoms"].symbols == atoms.symbols).all()

    sections = ["atoms", "eqpoints", "freqs", "kpoints", "mode_symmetries", "modes"]

    for key in sections:
        assert key in ph_results["results"][(0, 0, 0)]


def test_ph_job_list_to_do(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = bulk("Li")

    pp_dir = Path(__file__).parent

    copy_decompress_files([pp_dir / "Li.upf"], tmp_path)

    SETTINGS.ESPRESSO_PSEUDO = Path(__file__).parent

    input_data = {
        "calculation": "scf",
        "occupations": "smearing",
        "smearing": "cold",
        "degauss": 0.02,
        "mixing_mode": "TF",
        "mixing_beta": 0.7,
        "conv_thr": 1.0e-6,
    }

    ph_loose = {"tr2_ph": 1e-8, "qplot": True, "nat_todo": 1, "ldisp": True}

    pseudopotentials = {"Li": "Li.upf"}

    pw_results = static_job(
        atoms, input_data=input_data, pseudopotentials=pseudopotentials, kspacing=0.25
    )

    qpts = [(0, 0, 0, 1), (1 / 3, 0, 0, 1), (1 / 2, 0, 0, 1)]

    nat_todo = [1]

    ph_results = ph_job(
        input_data=ph_loose,
        copy_files=pw_results["dir_name"],
        qpts=qpts,
        nat_todo=nat_todo,
    )

    assert (0, 0, 0) in ph_results["results"]
    assert np.allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].positions,
        atoms.positions,
        atol=1.0e-4,
    )
    # ph.x cell param are not defined to a very high level of accuracy,
    # atol = 1.0e-3 is needed here...
    assert np.allclose(
        ph_results["results"][(0, 0, 0)]["atoms"].cell, atoms.cell, atol=1.0e-3
    )
    assert (ph_results["results"][(0, 0, 0)]["atoms"].symbols == atoms.symbols).all()

    sections = ["atoms", "eqpoints", "freqs", "kpoints", "mode_symmetries", "modes"]

    for key in sections:
        assert key in ph_results["results"][(0, 0, 0)]
