import os

import pytest
from ase.build import bulk
from jobflow import JobStore, run_locally
from maggma.stores import MemoryStore

from quacc.recipes.emt.core import relax_job, static_job
from quacc.recipes.emt.slabs import BulkToSlabsFlow

try:
    import jobflow as jf
except ImportError:
    jf = None

store = JobStore(MemoryStore())


@pytest.mark.skipif(jf is None, reason="This test requires jobflow")
def test_emt():
    atoms = bulk("Cu")

    job = jf.job(static_job)(atoms)
    run_locally(job, store=store, ensure_success=True)

    job = jf.job(relax_job)(atoms)
    run_locally(job, store=store, ensure_success=True)
