import pytest
from maggma.stores import MemoryStore

from quacc import SETTINGS, flow, job, subflow

jf = pytest.importorskip("jobflow")

pytestmark = pytest.mark.skipif(
    SETTINGS.WORKFLOW_ENGINE != "jobflow",
    reason="Jobflow needs to be the workflow engine",
)

STORE = jf.JobStore(MemoryStore())


def test_jobflow_decorators(tmpdir):
    tmpdir.chdir()

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    assert not isinstance(add, jf.Job)
    assert not isinstance(mult, jf.Job)
    assert hasattr(add, "original")
    assert hasattr(mult, "original")
    assert isinstance(add(1, 2), jf.Job)
    assert isinstance(mult(1, 2), jf.Job)


def test_jobflow_decorators_args(tmpdir):
    tmpdir.chdir()

    @job()
    def add(a, b):
        return a + b

    @job()
    def mult(a, b):
        return a * b

    assert not isinstance(add, jf.Job)
    assert not isinstance(mult, jf.Job)
    assert hasattr(add, "original")
    assert hasattr(mult, "original")
    assert isinstance(add(1, 2), jf.Job)
    assert isinstance(mult(1, 2), jf.Job)
