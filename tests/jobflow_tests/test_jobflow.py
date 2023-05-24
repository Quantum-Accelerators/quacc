import os
from shutil import rmtree

import jobflow as jf
import pytest
from ase.build import bulk
from maggma.stores import MemoryStore

from quacc.recipes.emt.core import relax_job, static_job


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
        if "quacc-tmp" in f or f == "tmp_dir":
            rmtree(f)


@pytest.mark.skipif(jf is None, reason="This test requires jobflow")
def test_emt():
    store = jf.JobStore(MemoryStore())

    atoms = bulk("Cu")

    # Test inddividual jobs
    job = jf.job(static_job)(atoms)
    jf.run_locally(job, store=store, ensure_success=True)

    job = jf.job(relax_job)(atoms)
    jf.run_locally(job, store=store, ensure_success=True)

    # Test a Flow
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(static_job)(job1.output["atoms"])

    workflow = jf.Flow([job1, job2])
    jf.run_locally(workflow, create_folders=True, ensure_success=True)
