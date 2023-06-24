import jobflow as jf
import pytest
from ase.build import bulk

from quacc.recipes.emt.core import relax_job, static_job

try:
    import fireworks
except ImportError:
    fireworks = None


@pytest.mark.skipif(fireworks is None, reason="This test requires fireworks")
def test_fireworks():
    from jobflow.managers.fireworks import flow_to_workflow, job_to_firework

    atoms = bulk("Cu")

    # Test fireworks creation
    job1 = jf.job(relax_job)(atoms)
    job2 = jf.job(static_job)(job1.output["atoms"])

    workflow = jf.Flow([job1, job2])
    job_to_firework(job1)
    flow_to_workflow(workflow)
