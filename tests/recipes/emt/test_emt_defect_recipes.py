import os
from shutil import rmtree

import pytest
from ase.build import bulk

try:
    from quacc.recipes.emt.defects import bulk_to_defects_flow
except ImportError:
    bulk_to_defects_flow = None


def teardown_module():
    for f in os.listdir(os.getcwd()):
        if (
            f.endswith(".log")
            or f.endswith(".pckl")
            or f.endswith(".traj")
            or f.endswith(".out")
            or ".gz" in f
        ):
            os.remove(f)
        if "quacc-tmp" in f or "job_" in f or f == "tmp_dir":
            if os.path.islink(f):
                os.unlink(f)
            else:
                rmtree(f)


@pytest.mark.skipif(bulk_to_defects_flow is None, reason="requires quacc[defects]")
def test_bulk_to_defects_flow():
    atoms = bulk("Cu")
    output = bulk_to_defects_flow(atoms, defect_relax_kwargs={"opt_swaps": {"fmax": 5}})
    assert len(output) == 1
    assert len(output[0]["atoms"]) == 107
    # TODO: Add a lot of specific assertions to test that the output is correct
    # This will involve writing multiple tests for various combinations of inputs.
    # Refer to the `test_emt_recipes.py` file for examples. In general, the more
    # tests you have, the more robust your code will be to changes in either quacc
    # itself or in any dependencies quacc relies on, e.g. pymatgen.
    # In general, it is also good to ensure that as many lines of your added code
    # are covered by tests as possible.
