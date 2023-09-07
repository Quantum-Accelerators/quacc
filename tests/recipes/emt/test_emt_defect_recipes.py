import os

import pytest
from ase.build import bulk

try:
    from quacc.recipes.emt.defects import bulk_to_defects_flow
except ImportError:
    bulk_to_defects_flow = None


@pytest.mark.skipif(bulk_to_defects_flow is None, reason="requires quacc[defects]")
def test_bulk_to_defects_flow(tmpdir):
    os.chdir(tmpdir)

    atoms = bulk("Cu")
    output = bulk_to_defects_flow(atoms, defect_relax_kwargs={"opt_swaps": {"fmax": 5}})
    assert len(output) == 1
    assert len(output[0]["atoms"]) == 107

    os.chdir(tmpdir)

    atoms = bulk("Cu")
    output = bulk_to_defects_flow(
        atoms, run_static=None, defect_relax_kwargs={"opt_swaps": {"fmax": 5}}
    )
    assert len(output) == 1
    assert len(output[0]["atoms"]) == 107
