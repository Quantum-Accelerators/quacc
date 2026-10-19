from __future__ import annotations

import os
from pathlib import Path

import pytest
from ase.build import bulk, molecule
from ase.calculators.gulp import GULP

from quacc import change_settings
from quacc.recipes.gulp._base import run_and_summarize
from quacc.recipes.gulp.core import relax_job, static_job


def test_static_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")

    output = static_job(atoms)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert output["parameters"]["keywords"] == "gfnff"
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = static_job(atoms, keywords={"gwolf": True})
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = static_job(atoms, use_gfnff=False)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    atoms = bulk("Cu") * (2, 2, 2)
    output = static_job(atoms)
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = static_job(atoms, keywords={"#gwolf"})
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = static_job(atoms, use_gfnff=False)
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]


def test_relax_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")

    output = relax_job(atoms)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = relax_job(atoms, keywords={"gwolf": True})
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = relax_job(atoms, relax_cell=True, use_gfnff=False)
    assert output["molecule_metadata"]["natoms"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    atoms = bulk("Cu") * (2, 2, 2)
    output = relax_job(atoms, relax_cell=True)
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" in output["parameters"]["keywords"]
    assert "conv" not in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = relax_job(atoms, keywords={"gwolf": True})
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = relax_job(atoms, relax_cell=True, use_gfnff=False)
    assert output["structure_metadata"]["nsites"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" in output["parameters"]["keywords"]
    assert "conv" not in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]


def test_envvars(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    atoms = molecule("H2O")

    with change_settings({"GULP_LIB": str(Path("/path/to/lib"))}):
        assert static_job(atoms)
        assert os.environ.get("GULP_LIB") == str(Path("/path/to/lib"))


def test_unconverged_optimization(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(GULP, "get_opt_state", lambda self: False)

    with (
        change_settings({"CHECK_CONVERGENCE": True}),
        pytest.raises(RuntimeError, match="Optimization did not converge"),
    ):
        run_and_summarize(molecule("H2O"), keyword_defaults=["opti"])
