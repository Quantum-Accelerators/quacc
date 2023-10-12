import pytest

prefect = pytest.importorskip("prefect")


@pytest.fixture(autouse=True, scope="session")
def prefect_test_fixture():
    from prefect.testing.utilities import prefect_test_harness

    with prefect_test_harness():
        yield


def test_prefect_decorators(tmpdir):
    tmpdir.chdir()
    from quacc import flow, job, subflow

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

    @flow
    def add_flow(a, b):
        return add(a, b)

    assert add_flow(1, 2).result() == 3
    assert workflow(1, 2, 3).result() == 9
    results = dynamic_workflow(1, 2, 3)
    assert [result.result() for result in results] == [6, 6, 6]
