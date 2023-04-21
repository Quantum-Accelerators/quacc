import os
from pathlib import Path
from shutil import copy, rmtree

from ase.build import bulk, molecule
import covalent as ct
from quacc.recipes.gulp.core import RelaxJob, StaticJob

FILE_DIR = Path(__file__).resolve().parent
GULP_DIR = os.path.join(FILE_DIR, "gulp_run")


def setup_module():
    for f in os.listdir(GULP_DIR):
        copy(os.path.join(GULP_DIR, f), os.path.join(os.getcwd(), f))


def teardown_module():
    for f in os.listdir(GULP_DIR):
        if os.path.exists(os.path.join(os.getcwd(), f)):
            os.remove(os.path.join(os.getcwd(), f))
    for f in os.listdir(os.getcwd()):
        if "quacc-tmp" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


def test_static_Job():
    atoms = molecule("H2O")

    output = StaticJob(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["keywords"] == "gfnff"
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = StaticJob(atoms, keyword_swaps={"gwolf": True})
    assert output["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = StaticJob(atoms, gfnff=False)
    assert output["natoms"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    atoms = bulk("Cu") * (2, 2, 2)
    output = StaticJob(atoms)
    assert output["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = StaticJob(atoms, keyword_swaps={"gwolf": False})
    assert output["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = StaticJob(atoms, gfnff=False)
    assert output["natoms"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]


def test_relax_Job():
    atoms = molecule("H2O")

    output = RelaxJob(atoms)
    assert output["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = RelaxJob(atoms, volume_relax=False, keyword_swaps={"gwolf": True})
    assert output["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = RelaxJob(atoms, gfnff=False)
    assert output["natoms"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    atoms = bulk("Cu") * (2, 2, 2)
    output = RelaxJob(atoms)
    assert output["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" in output["parameters"]["keywords"]
    assert "conv" not in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = RelaxJob(atoms, volume_relax=False, keyword_swaps={"gwolf": True})
    assert output["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = RelaxJob(atoms, gfnff=False)
    assert output["natoms"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" in output["parameters"]["keywords"]
    assert "conv" not in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]
