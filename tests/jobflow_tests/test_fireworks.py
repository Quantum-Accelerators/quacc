import pytest
from ase.build import bulk

from quacc import SETTTINGS
from quacc.recipes.emt.core import relax_job, static_job

try:
    import fireworks
    import jobflow as jf

    jf = jf if SETTINGS.WORKFLOW_ENGINE == "jobflow" else None

except ImportError:
    fireworks = None


@pytest.mark.skipif(fireworks is None, reason="This test requires fireworks")
def test_fireworks():
    from jobflow.managers.fireworks import flow_to_workflow, job_to_firework

    atoms = bulk("Cu")

    # Test fireworks creation
    job1 = relax_job(atoms)
    job2 = static_job(job1.output)

    workflow = jf.Flow([job1, job2])
    assert job_to_firework(job1)
    assert flow_to_workflow(workflow)
