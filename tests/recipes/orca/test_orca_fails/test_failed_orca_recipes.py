import pytest
from ase.build import molecule


def test_static_job(tmpdir):
    from quacc.recipes.orca.core import relax_job

    tmpdir.chdir()

    atoms = molecule("H2")
    with pytest.raises(ValueError):
        relax_job(atoms, 0, 1)
