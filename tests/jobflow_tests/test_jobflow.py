import os
from shutil import rmtree

import pytest
from ase.build import bulk

from quacc.recipes.emt.core import relax_job, static_job

try:
    import fireworks
    import jobflow as jf
    import maggma

    missing_imports = False
except ImportError:
    missing_imports = True


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


@pytest.mark.skipif(
    missing_imports, reason="This test requires jobflow, fireworks, and maggma"
)
def test_emt():
    import jobflow as jf
    from jobflow.managers.fireworks import flow_to_workflow, job_to_firework
    from maggma.stores import MemoryStore

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

    # Test fireworks creation
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(static_job)(job1.output["atoms"])

    workflow = jf.Flow([job1, job2])
    job_to_firework(job1)
    flow_to_workflow(workflow)
