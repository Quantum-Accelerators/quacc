from __future__ import annotations

import pytest

pytest.importorskip("aiida_workgraph")


from quacc import flow, job


def test_aiida_decorators():
    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    assert workflow.run(1, 2, 3) == 9
