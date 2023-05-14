import os

import pytest
from ase.build import bulk
from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.slabs import BulkToSlabsFlow

try:
    import covalent as ct
except ImportError:
    ct = None


@pytest.mark.skipif(
    ct is None,
    reason="covalent must be installed. Try pip install covalent",
)
def test_covalent_config():
    ct_config = ct.get_config()
    if ct_config["executors"].get("slurm", None):
        assert ct_config["executors"]["slurm"].get("use_srun", "true") == "false"
        assert (
            ct_config["executors"]["slurm"].get("create_unique_workdir", "false")
            == "true"
        )


@pytest.mark.skipif(
    ct is None,
    reason="covalent must be installed. Try pip install covalent",
)
def test_emt():
    atoms = bulk("Cu")
    workflow = ct.lattice(ct.electron(static_job))
    dispatch_id = ct.dispatch(workflow)(atoms)
    result = ct.get_result(dispatch_id, wait=True)
    assert result.status == "COMPLETED"
