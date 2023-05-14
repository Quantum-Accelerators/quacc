import os

import pytest
from ase.build import bulk
from jobflow.managers.local import run_locally

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.slabs import BulkToSlabsFlow

try:
    import jobflow as jf
except ImportError:
    jf = None


@pytest.mark.skipif(
    jf is None or os.environ.get("GITHUB_ACTIONS", False) is False,
    reason="This test is only meant to be run on GitHub Actions",
)
def test_emt():
    atoms = bulk("Cu")

    job = jf.job(static_job)(atoms)
    run_locally(job, ensure_success=True)

    job = jf.job(relax_job)(atoms)
    run_locally(job, ensure_success=True)
