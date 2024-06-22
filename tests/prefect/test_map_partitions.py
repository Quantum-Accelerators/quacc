from __future__ import annotations

import numpy as np
import pytest

from quacc import job, flow
from quacc.wflow_tools.job_patterns import (
    kwarg_map,
    map_partition,
    map_partitioned_lists,
    partition,
)


def test_partition():
    @job
    def simple_list():
        return list(range(10))

    @flow
    def test_flow():
        return partition(simple_list(), 3)

    partitioned_list = test_flow().result()
    assert len(partitioned_list) == 3
    np.testing.assert_allclose(partitioned_list[0], [0, 1, 2, 3])


def test_kwarg_map():
    def test_fun(a, b):
        return {"a": a, "b": b}

    assert kwarg_map(test_fun, a=[1, 2], b=["c", "d"]) == [
        {"a": 1, "b": "c"},
        {"a": 2, "b": "d"},
    ]
    assert kwarg_map(test_fun, unmapped_kwargs={"b": 1}, a=[1, 2]) == [
        {"a": 1, "b": 1},
        {"a": 2, "b": 1},
    ]
    with pytest.raises(AssertionError):
        kwarg_map(test_fun, a=[1, 2, 3], b=[1, 2])


def test_map_partitioned_lists():
    @job
    def simple_list():
        return list(range(10))

    def testfun(a, const=2):
        return a * const

    @flow
    def test_flow(func, num_partitions, unmapped_kwargs=None):
        return map_partitioned_lists(
            func,
            unmapped_kwargs=unmapped_kwargs,
            a=partition(simple_list(), num_partitions),
            num_partitions=num_partitions,
        )

    num_partitions = 4
    assert test_flow(testfun, num_partitions)[0].result() == [0, 2, 4]
    assert test_flow(testfun, num_partitions, unmapped_kwargs={"const": 3})[
        0
    ].result() == [0, 3, 6]
    num_partitions = 2
    assert test_flow(testfun, num_partitions)[0].result() == [0, 2, 4, 6, 8]
    assert test_flow(testfun, num_partitions, unmapped_kwargs={"const": 3})[
        1
    ].result() == [15, 18, 21, 24, 27]


def test_map_partition():
    def test_fun(a, b):
        return {"a": a, "b": b}

    @flow
    def test_flow1():
        return map_partition(test_fun, a=[1, 2], b=["c", "d"])

    assert test_flow1().result() == [{"a": 1, "b": "c"}, {"a": 2, "b": "d"}]

    @flow
    def test_flow2():
        return map_partition(test_fun, unmapped_kwargs={"b": 1}, a=[1, 2])

    assert test_flow2().result() == [{"a": 1, "b": 1}, {"a": 2, "b": 1}]
