from __future__ import annotations

import itertools
from importlib.util import find_spec
from typing import TYPE_CHECKING

from quacc.wflow_tools.customizers import strip_decorator
from quacc.wflow_tools.decorators import job

has_fairchem = bool(find_spec("fairchem"))

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from ase.atoms import Atoms


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


def map_partitioned_lists_fairchembatch(
    func: Callable,
    num_partitions: int,
    name_or_path: str,
    task_name: str | None = None,
    inference_settings: str = "default",
    device: str | None = None,
    unmapped_kwargs: dict[str, Any] | None = None,
    batcher_kwargs: dict[str, Any] | None = None,
    **mapped_kwargs: list[list[Any]],
) -> list[Any]:
    """
    Given list-of-lists parameters (partitioned data), apply func to each element
    using FAIRChem batched inference.

    This is the FAIRChem batched equivalent of `map_partitioned_lists`. Each partition
    is dispatched as a separate job via `map_partition_fairchembatch`, which uses
    Ray Serve for efficient GPU batching within each partition.

    Parameters
    ----------
    func
        The function to map. The function will receive `predict_unit` and `task_name`
        automatically injected into its kwargs (along with `method="fairchem"`).
    num_partitions
        The number of partitions (length of each list in mapped_kwargs).
    name_or_path
        A model name from fairchem.core.pretrained.available_models or a path
        to the checkpoint file (e.g., "uma-s-1", "uma-m-1").
    task_name
        Task name (e.g., 'omat', 'omol', 'oc20', 'odac', 'omc').
    inference_settings
        Settings for inference. Can be "default" or "turbo".
    device
        Optional torch device to load the model onto (e.g., 'cuda', 'cpu').
    unmapped_kwargs
        Dictionary of kwargs to pass to func that shouldn't be mapped.
    batcher_kwargs
        Additional kwargs to pass to get_inference_batcher.
    **mapped_kwargs
        kwargs of the form key=list[list[...]] that should be mapped over.

    Returns
    -------
    list[Any]
        List of results (one per partition), each containing results for that partition.

    Examples
    --------
    >>> from quacc.recipes.mlp.core import relax_job
    >>> from quacc.wflow_tools.job_patterns import (
    ...     map_partitioned_lists_fairchembatch,
    ...     partition,
    ...     unpartition,
    ... )
    >>>
    >>> num_partitions = 4
    >>> partitioned_atoms = partition(my_atoms_list, num_partitions)
    >>>
    >>> results_partitioned = map_partitioned_lists_fairchembatch(
    ...     relax_job,
    ...     num_partitions,
    ...     name_or_path="uma-s-1",
    ...     task_name="omat",
    ...     atoms=partitioned_atoms,
    ... )
    >>>
    >>> all_results = unpartition(results_partitioned)
    """
    return [
        map_partition_fairchembatch(
            strip_decorator(func),
            name_or_path=name_or_path,
            task_name=task_name,
            inference_settings=inference_settings,
            device=device,
            unmapped_kwargs=unmapped_kwargs,
            batcher_kwargs=batcher_kwargs,
            **{k: mapped_kwargs[k][i] for k in mapped_kwargs},
        )
        for i in range(num_partitions)
    ]


@job
def map_partition_fairchembatch(
    func: Callable,
    name_or_path: str,
    task_name: str | None = None,
    inference_settings: str = "default",
    device: str | None = None,
    unmapped_kwargs: dict[str, Any] | None = None,
    batcher_kwargs: dict[str, Any] | None = None,
    **mapped_kwargs: list[Any],
) -> list[Any]:
    """
    Job to apply a function to each set of kwargs using FAIRChem batched inference.

    This is the FAIRChem batched equivalent of `map_partition`. It processes a single
    partition of inputs, submitting calculations concurrently via threads while
    the underlying Ray Serve batches requests to the GPU.

    The InferenceBatcher is cached based on the model configuration, so repeated
    calls with the same model will reuse the same Ray Serve deployment.

    Parameters
    ----------
    func
        The function to map. The function will receive `predict_unit` and `task_name`
        automatically injected into its kwargs (along with `method="fairchem"`).
    name_or_path
        A model name from fairchem.core.pretrained.available_models or a path
        to the checkpoint file (e.g., "uma-s-1", "uma-m-1").
    task_name
        Task name (e.g., 'omat', 'omol', 'oc20', 'odac', 'omc').
    inference_settings
        Settings for inference. Can be "default" (general purpose) or "turbo"
        (optimized for speed but requires fixed atomic composition).
    device
        Optional torch device to load the model onto (e.g., 'cuda', 'cpu').
    unmapped_kwargs
        Dictionary of kwargs to pass to func that shouldn't be mapped.
        These are merged with the FAIRChem-specific kwargs (predict_unit, task_name).
    batcher_kwargs
        Additional kwargs to pass to get_inference_batcher. Available options:
        - max_batch_size: Maximum number of atoms in a batch
        - batch_wait_timeout_s: Max time to wait for batch
        - num_replicas: Number of Ray Serve replicas
        - concurrency_backend_options: Dict with options like {"max_workers": N}
    **mapped_kwargs
        kwargs of the form key=list[...] that should be mapped over.
        All lists must have the same length.

    Returns
    -------
    list[Any]
        List of results from calling func for each set of mapped kwargs.
    """
    if not has_fairchem:
        raise ImportError(
            "fairchem must be installed to use map_partition_fairchembatch. "
            "Run pip install fairchem-core."
        )

    from quacc.recipes.mlp._base import get_inference_batcher

    # Validate mapped_kwargs lengths
    if not mapped_kwargs:
        raise ValueError("At least one mapped kwarg must be provided")

    all_lens = [len(v) for v in mapped_kwargs.values()]
    n_elements = all_lens[0]
    if not all(n_elements == le for le in all_lens):
        raise AssertionError(f"Inconsistent lengths: {all_lens}")

    # Get or create cached batcher
    batcher_kwargs = batcher_kwargs or {}
    batcher = get_inference_batcher(
        name_or_path=name_or_path,
        task_name=task_name,
        inference_settings=inference_settings,
        device=device,
        **batcher_kwargs,
    )

    # Prepare base kwargs for all calls
    unmapped_kwargs = unmapped_kwargs or {}
    base_kwargs = {
        "method": "fairchem",
        "predict_unit": batcher.batch_predict_unit,
        "task_name": task_name,
        **unmapped_kwargs,
    }

    # Strip decorator from func if needed
    raw_func = strip_decorator(func)

    # Submit all tasks concurrently using the batcher's thread pool
    futures = []
    for i in range(n_elements):
        # Build kwargs for this specific call
        call_kwargs = {**base_kwargs}
        for key, values in mapped_kwargs.items():
            call_kwargs[key] = values[i]

        future = batcher.executor.submit(raw_func, **call_kwargs)
        futures.append(future)

    # Wait for all results
    return [future.result() for future in futures]
