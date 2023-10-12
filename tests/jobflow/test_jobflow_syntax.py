import pytest

jf = pytest.importorskip("jobflow")

from quacc import SETTINGS

DEFAULT_SETTINGS = SETTINGS.copy()


def setup_module():
    SETTINGS.WORKFLOW_ENGINE = "jobflow"


def teardown_module():
    from quacc import SETTINGS

    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE


def test_jobflow_decorators(tmpdir):
    from quacc import flow, job, subflow

    tmpdir.chdir()

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    assert not isinstance(add, jf.Job)
    assert not isinstance(mult, jf.Job)
    assert hasattr(add, "original")
    assert hasattr(mult, "original")
    assert isinstance(add(1, 2), jf.Job)
    assert isinstance(mult(1, 2), jf.Job)
    assert isinstance(workflow(1, 2, 3), jf.Job)
    assert isinstance(add_distributed([1, 2, 3], 4)[0], jf.Job)


def test_jobflow_decorators_args(tmpdir):
    from quacc import flow, job, subflow

    tmpdir.chdir()

    @job()
    def add(a, b):
        return a + b

    @job()
    def mult(a, b):
        return a * b

    @subflow()
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow()
    def workflow(a, b, c):
        return mult(add(a, b), c)

    assert not isinstance(add, jf.Job)
    assert not isinstance(mult, jf.Job)
    assert hasattr(add, "original")
    assert hasattr(mult, "original")
    assert isinstance(add(1, 2), jf.Job)
    assert isinstance(mult(1, 2), jf.Job)
    assert isinstance(workflow(1, 2, 3), jf.Job)
    assert isinstance(add_distributed([1, 2, 3], 4)[0], jf.Job)
