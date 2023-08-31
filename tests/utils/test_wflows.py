import os

import pytest

from quacc import SETTINGS, flow, job, subflow

try:
    import covalent as ct
except ImportError:
    ct = None
try:
    import parsl
except ImportError:
    parsl = None
try:
    import jobflow
except ImportError:
    jobflow = None

DEFAULT_SETTINGS = SETTINGS.copy()


def teardown_module():
    SETTINGS.WORKFLOW_ENGINE = DEFAULT_SETTINGS.WORKFLOW_ENGINE


def test_decorators(tmpdir):
    tmpdir.chdir()

    SETTINGS.WORKFLOW_ENGINE = "local"

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2) == 3
    assert not hasattr(add, "electron_object")
    assert mult(1, 2) == 2
    assert not hasattr(mult, "electron_object")
    assert workflow(1, 2, 3) == 9
    assert not hasattr(workflow, "electron_object")
    assert dynamic_workflow(1, 2, 3) == [6, 6, 6]
    assert not hasattr(dynamic_workflow, "electron_object")
    assert add_distributed([1, 2, 3], 4) == [5, 6, 7]
    assert not hasattr(add_distributed, "electron_object")


def test_decorators_args(tmpdir):
    tmpdir.chdir()

    SETTINGS.WORKFLOW_ENGINE = "local"

    @job()
    def add(a, b):
        return a + b

    @job()
    def mult(a, b):
        return a * b

    @job()
    def make_more(val):
        return [val] * 3

    @subflow()
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow()
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow()
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2) == 3
    assert not hasattr(add, "electron_object")
    assert mult(1, 2) == 2
    assert not hasattr(mult, "electron_object")
    assert workflow(1, 2, 3) == 9
    assert not hasattr(workflow, "electron_object")
    assert dynamic_workflow(1, 2, 3) == [6, 6, 6]
    assert not hasattr(dynamic_workflow, "electron_object")
    assert add_distributed([1, 2, 3], 4) == [5, 6, 7]
    assert not hasattr(add_distributed, "electron_object")


@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS", False) is False or not ct,
    reason="This test requires Covalent and to be run on GitHub",
)
def test_covalent_decorators(tmpdir):
    tmpdir.chdir()

    SETTINGS.WORKFLOW_ENGINE = "covalent"

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2) == 3
    assert mult(1, 2) == 2
    assert ct.get_result(workflow(1, 2, 3), wait=True).result == 9
    assert ct.get_result(dynamic_workflow(1, 2, 3), wait=True).result == [6, 6, 6]


@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS", False) is False or not ct,
    reason="This test requires Covalent and to be run on GitHub",
)
def test_covalent_decorators_args(tmpdir):
    tmpdir.chdir()

    SETTINGS.WORKFLOW_ENGINE = "covalent"

    @job(executor="local")
    def add(a, b):
        return a + b

    @job(executor="local")
    def mult(a, b):
        return a * b

    @job(executor="local")
    def make_more(val):
        return [val] * 3

    @subflow(executor="local")
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow(executor="local")
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow(executor="local")
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2) == 3
    assert mult(1, 2) == 2
    assert ct.get_result(workflow(1, 2, 3), wait=True).result == 9
    assert ct.get_result(dynamic_workflow(1, 2, 3), wait=True).result == [6, 6, 6]


@pytest.mark.skipif(parsl is None, reason="Parsl not installed")
def test_parsl_decorators(tmpdir):
    tmpdir.chdir()
    SETTINGS.WORKFLOW_ENGINE = "parsl"

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2).result() == 3
    assert mult(1, 2).result() == 2
    assert workflow(1, 2, 3).result() == 9
    assert dynamic_workflow(1, 2, 3).result() == [6, 6, 6]


@pytest.mark.skipif(parsl is None, reason="Parsl not installed")
def test_parsl_decorators_args(tmpdir):
    tmpdir.chdir()
    SETTINGS.WORKFLOW_ENGINE = "parsl"

    @job()
    def add(a, b):
        return a + b

    @job()
    def mult(a, b):
        return a * b

    @job()
    def make_more(val):
        return [val] * 3

    @subflow()
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow()
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow()
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2).result() == 3
    assert mult(1, 2).result() == 2
    assert workflow(1, 2, 3).result() == 9
    assert dynamic_workflow(1, 2, 3).result() == [6, 6, 6]


@pytest.mark.skipif(jobflow is None, reason="Jobflow not installed")
def test_jobflow_decorators(tmpdir):
    tmpdir.chdir()

    SETTINGS.WORKFLOW_ENGINE = "jobflow"
    from jobflow import Job

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

    assert not isinstance(add, Job)
    assert not isinstance(mult, Job)
    assert isinstance(add(1, 2), Job)
    assert isinstance(mult(1, 2), Job)
    assert isinstance(workflow(1, 2, 3), Job)
    assert isinstance(add_distributed([1, 2, 3], 4)[0], Job)


@pytest.mark.skipif(jobflow is None, reason="Jobflow not installed")
def test_jobflow_decorators_args(tmpdir):
    tmpdir.chdir()

    SETTINGS.WORKFLOW_ENGINE = "jobflow"
    from jobflow import Job

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

    assert not isinstance(add, Job)
    assert not isinstance(mult, Job)
    assert isinstance(add(1, 2), Job)
    assert isinstance(mult(1, 2), Job)
    assert isinstance(workflow(1, 2, 3), Job)
    assert isinstance(add_distributed([1, 2, 3], 4)[0], Job)
