import pytest

from quacc import flow, job, subflow

try:
    import prefect

except ImportError:
    prefect = None


@pytest.mark.skipif(
    prefect is None,
    reason="This test requires Prefect",
)
def test_prefect_decorators(tmpdir):
    tmpdir.chdir()

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
