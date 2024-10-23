from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

from quacc.wflow_tools.customizers import strip_decorator
from quacc.wflow_tools.decorators import job

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any


@job
def partition(list_to_partition: list, num_partitions: int) -> list[Any]:
    """
    Given a list, partition it into n roughly equal lists

    Parameters
    ----------
    list_to_partition
        the list to partition
    num_partitions
        the number of partitions to output

    Returns
    -------
    list[Any]
        n lists constructed from a
    """
    k, m = divmod(len(list_to_partition), num_partitions)
    return [
        list_to_partition[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)]
        for i in range(num_partitions)
    ]


def map_partitioned_lists(
    func: Callable,
    num_partitions: int,
    unmapped_kwargs: dict[str, Any] | None = None,
    **mapped_kwargs: list[list[Any]],
) -> list[Any]:
    """
    Given list-of-lists parameters (say a list of batches that we want to map over),
    apply func to each element of each list

    For example:

    ```python
    @job
    def testjob(**kwargs):
        print(kwargs)


    @flow
    def testflow():
        num_partitions = 2
        result = map_partitioned_lists(
            testjob,
            num_partitions,
            test_arg_1=partition([1, 2, 3, 4, 5], num_partitions),
            test_arg_2=partition(["a", "b", "c", "d", "e"], num_partitions),
        )


    testflow()
    ```

    should yield

    ```python
    {"test_arg_1": 1, "test_arg_2": "a"}
    {"test_arg_1": 2, "test_arg_2": "b"}
    {"test_arg_1": 3, "test_arg_2": "c"}
    {"test_arg_1": 4, "test_arg_2": "d"}
    {"test_arg_1": 5, "test_arg_2": "e"}
    ```

    regardless of the number of partitions.


    Parameters
    ----------
    func
        The function to map.
    num_partitions
        the length of each kwarg in mapped_kwargs
    unmapped_kwargs
        Dictionary of kwargs to pass to func that shouldn't be mapped
    mapped_kwargs
        kwargs of the form key=list[...] that should be mapped over

    Returns
    -------
    list[Any]
        list of results from calling func(**(mapped_kwargs | unmapped_kwargs)) for each
        kwargs in mapped_kwargs
    """

    return [
        map_partition(
            strip_decorator(func),
            unmapped_kwargs=unmapped_kwargs,
            **{k: mapped_kwargs[k][i] for k in mapped_kwargs},
        )
        for i in range(num_partitions)
    ]


@job
def map_partition(
    func: Callable, unmapped_kwargs: dict[str, Any] | None = None, **mapped_kwargs
) -> list[Any]:
    """
    Job to apply a function to each set of elements in mapped_kwargs.

    Parameters
    ----------
    func
        The function to map.
    unmapped_kwargs
        Dictionary of kwargs to pass to func that shouldn't be mapped
    mapped_kwargs
        kwargs of the form key=list[...] that should be mapped over

    Returns
    -------
    list[Any]
        list of results from calling func(**mapped_kwargs, **unmapped_kwargs) for each
        kwargs in mapped_kwargs
    """
    return kwarg_map(func, unmapped_kwargs=unmapped_kwargs, **mapped_kwargs)


def kwarg_map(
    func: Callable, unmapped_kwargs: dict[str, Any] | None = None, **mapped_kwargs
) -> list[Any]:
    """
    A helper function for when you want to construct a chain of objects with individual arguments for each one.  Can
    be easier to read than a list expansion.

    Adapted from https://stackoverflow.com/a/36575917 (CC-by-SA 3.0)

    Parameters
    ----------
    func
        The function to map.
    unmapped_kwargs
        Dictionary of kwargs to pass to func that shouldn't be mapped
    mapped_kwargs
        kwargs of the form key=list[...] that should be mapped over

    Returns
    -------
    list[Any]
        List of results from calling func(**mapped_kwargs, **unmapped_kwargs) for each
        kwargs in mapped_kwargs
    """
    unmapped_kwargs = unmapped_kwargs or {}

    all_lens = [len(v) for v in mapped_kwargs.values()]
    n_elements = all_lens[0]
    if not all(n_elements == le for le in all_lens):
        raise AssertionError(f"Inconsistent lengths: {all_lens}")
    return [
        func(**{k: v[i] for k, v in iter(mapped_kwargs.items())}, **unmapped_kwargs)
        for i in range(n_elements)
    ]


@job
def unpartition(lists_to_combine: list[list[Any]]) -> list[Any]:
    """
    Given a partitioned list (list of lists), recombine
    it to a single list

    Parameters
    ----------
    lists_to_combine
        the list of lists to recombine

    Returns
    -------
    list[Any]
        a single recombined list
    """
    return list(itertools.chain(*lists_to_combine))
