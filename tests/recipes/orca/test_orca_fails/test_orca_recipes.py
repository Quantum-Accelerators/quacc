import pytest
from ase.build import molecule

from quacc.recipes.orca.core import relax_job


def test_static_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2")
    with pytest.raises(ValueError):
        relax_job(atoms)
