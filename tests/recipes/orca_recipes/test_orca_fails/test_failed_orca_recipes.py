import pytest
from ase.build import molecule

from quacc import SETTINGS

pytestmark = pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE != "local",
    reason="Need to use local as workflow manager to run this test.",
)


def test_static_job(tmpdir):
    from quacc.recipes.orca.core import relax_job

    tmpdir.chdir()

    atoms = molecule("H2")
    with pytest.raises(ValueError):
        relax_job(atoms, 0, 1)
