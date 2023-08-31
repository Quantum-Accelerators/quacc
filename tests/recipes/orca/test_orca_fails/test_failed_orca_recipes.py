import pytest
from ase.build import molecule

from quacc import SETTINGS
from quacc.recipes.orca.core import relax_job


@pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE not in {"local", "covalent"},
    reason="This test suite is for regular function execution only",
)
def test_static_job(tmpdir):
    tmpdir.chdir()

    atoms = molecule("H2")
    with pytest.raises(ValueError):
        relax_job(atoms)
