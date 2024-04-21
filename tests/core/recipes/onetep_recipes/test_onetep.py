from __future__ import annotations

from pathlib import Path

from ase.build import molecule
from ase.optimize import LBFGS

from quacc.recipes.onetep.core import ase_relax_job, static_job
from quacc.utils.files import copy_decompress_files

DATA_DIR = Path(__file__).parent / "data"


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files(DATA_DIR, ["H.usp.gz"], tmp_path)

    keywords = {"pseudo_path": str(tmp_path)}
    pseudopotentials = {"H": "H.usp"}

    atoms = molecule("H2")
    atoms.set_cell([10, 10, 10])
    atoms.pbc = True
    atoms.center()

    output = static_job(
        atoms,
        keywords=keywords,
        pseudopotentials=pseudopotentials,
        ngwf_radius=6.0,
        ngwf_count=1,
    )

    assert output["atoms"] == atoms
    assert output["parameters"]["keywords"] == {
        "output_detail": "verbose",
        "do_properties": True,
        "cutoff_energy": "600 eV",
        "task": "SinglePoint",
        "pseudo_path": str(tmp_path),
    }


def test_ase_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    copy_decompress_files(DATA_DIR, ["H.usp.gz"], tmp_path)

    keywords = {"pseudo_path": str(tmp_path)}
    pseudopotentials = {"H": "H.usp"}
    opt_params = {"fmax": 5.0, "max_steps": 2, "optimizer": LBFGS}
    atoms = molecule("H2")
    atoms.set_cell([10, 10, 10])
    atoms.pbc = True
    atoms.center()

    output = ase_relax_job(
        atoms,
        opt_params=opt_params,
        keywords=keywords,
        pseudopotentials=pseudopotentials,
        ngwf_radius=6.0,
        ngwf_count=1,
    )

    assert output["atoms"] == atoms
    assert output["parameters"]["keywords"] == {
        "output_detail": "verbose",
        "do_properties": True,
        "cutoff_energy": "600 eV",
        "task": "SinglePoint",
        "write_forces": True,
        "pseudo_path": str(tmp_path),
        "forces_output_detail": "verbose",
    }
    assert output["parameters_opt"]["max_steps"] == 2
    assert output["parameters_opt"]["optimizer"] == "LBFGS"
