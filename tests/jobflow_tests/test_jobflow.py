import os
from shutil import rmtree

import jobflow as jf
import pytest
from ase.build import bulk
from jobflow.managers.fireworks import flow_to_workflow, job_to_firework
from maggma.stores import MemoryStore

from quacc.recipes.emt.core import relax_job, static_job

try:
    import fireworks
except ImportError:
    fireworks = True


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
            rmtree(f)


def test_emt():
    store = jf.JobStore(MemoryStore())

    atoms = bulk("Cu")

    # Test individual jobs
    job = jf.job(static_job)(atoms)
    jf.run_locally(job, store=store, ensure_success=True)

    job = jf.job(relax_job)(atoms)
    jf.run_locally(job, store=store, ensure_success=True)

    # Test a Flow
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(static_job)(job1.output["atoms"])

    workflow = jf.Flow([job1, job2])
    jf.run_locally(workflow, create_folders=True, ensure_success=True)


@pytest.mark.skipif(fireworks is None, reason="This test requires fireworks")
def test_fireworks():
    atoms = bulk("Cu")

    # Test fireworks creation
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(static_job)(job1.output["atoms"])

    workflow = jf.Flow([job1, job2])
    job_to_firework(job1)
    flow_to_workflow(workflow)
