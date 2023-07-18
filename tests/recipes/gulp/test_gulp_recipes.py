from ase.build import bulk, molecule

from quacc.recipes.gulp.core import relax_job, static_job


def test_static_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")

    output = static_job(atoms)
    assert output["natoms"] == len(atoms)
    assert output["parameters"]["keywords"] == "gfnff"
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = static_job(atoms, keyword_swaps={"gwolf": True})
    assert output["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = static_job(atoms, use_gfnff=False)
    assert output["natoms"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    atoms = bulk("Cu") * (2, 2, 2)
    output = static_job(atoms)
    assert output["nsites"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = static_job(atoms, keyword_swaps={"gwolf": None})
    assert output["nsites"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = static_job(atoms, use_gfnff=False)
    assert output["nsites"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]


def test_relax_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2O")

    output = relax_job(atoms)
    assert output["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = relax_job(atoms, relax_cell=False, keyword_swaps={"gwolf": True})
    assert output["natoms"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" in output["parameters"]["options"]
    assert "output cif gulp.cif" not in output["parameters"]["options"]

    output = relax_job(atoms, use_gfnff=False)
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
    output = relax_job(atoms)
    assert output["nsites"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" in output["parameters"]["keywords"]
    assert "conv" not in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = relax_job(atoms, relax_cell=False, keyword_swaps={"gwolf": True})
    assert output["nsites"] == len(atoms)
    assert "gfnff" in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" not in output["parameters"]["keywords"]
    assert "conv" in output["parameters"]["keywords"]
    assert "gwolf" in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]

    output = relax_job(atoms, use_gfnff=False)
    assert output["nsites"] == len(atoms)
    assert "gfnff" not in output["parameters"]["keywords"]
    assert "opti" in output["parameters"]["keywords"]
    assert "conp" in output["parameters"]["keywords"]
    assert "conv" not in output["parameters"]["keywords"]
    assert "gwolf" not in output["parameters"]["keywords"]
    assert "dump every gulp.res" in output["parameters"]["options"]
    assert "output xyz gulp.xyz" not in output["parameters"]["options"]
    assert "output cif gulp.cif" in output["parameters"]["options"]
