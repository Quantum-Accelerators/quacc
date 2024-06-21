from __future__ import annotations

import numpy as np
import pytest

from quacc import job
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

    partitioned_list = partition(simple_list(), 3)
    assert len(partitioned_list) == 3
    assert np.allclose(partitioned_list[0], [0, 1, 2, 3])

    return


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
    return


def test_map_partitioned_lists():
    @job
    def simple_list():
        return list(range(10))

    def testfun(a, const=2):
        return a * const

    num_partitions = 4
    assert map_partitioned_lists(
        testfun,
        a=partition(simple_list(), num_partitions),
        num_partitions=num_partitions,
    )[0] == [0, 2, 4]
    assert map_partitioned_lists(
        testfun,
        unmapped_kwargs={"const": 3},
        a=partition(simple_list(), num_partitions),
        num_partitions=num_partitions,
    )[0] == [0, 3, 6]
    num_partitions = 2
    assert map_partitioned_lists(
        testfun,
        a=partition(simple_list(), num_partitions),
        num_partitions=num_partitions,
    )[0] == [0, 2, 4, 6, 8]
    assert map_partitioned_lists(
        testfun,
        unmapped_kwargs={"const": 3},
        a=partition(simple_list(), num_partitions),
        num_partitions=num_partitions,
    )[1] == [15, 18, 21, 24, 27]
    return


def test_map_partition():
    def test_fun(a, b):
        return {"a": a, "b": b}

    assert map_partition(test_fun, a=[1, 2], b=["c", "d"]) == [
        {"a": 1, "b": "c"},
        {"a": 2, "b": "d"},
    ]
    assert map_partition(test_fun, unmapped_kwargs={"b": 1}, a=[1, 2]) == [
        {"a": 1, "b": 1},
        {"a": 2, "b": 1},
    ]
    with pytest.raises(AssertionError):
        kwarg_map(test_fun, a=[1, 2, 3], b=[1, 2])
    return
