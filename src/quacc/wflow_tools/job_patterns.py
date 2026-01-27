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


@job
def map_partitioned_lists_fairchembatch(
    func: Callable,
    atoms_list: list[Atoms],
    fairchem_model: str,
    task_name: str | None = None,
    inference_settings: str = "default",
    device: str | None = None,
    unmapped_kwargs: dict[str, Any] | None = None,
    batcher_kwargs: dict[str, Any] | None = None,
    **mapped_kwargs: list[Any],
) -> list[Any]:
    """
    Apply a function to each atoms object using FAIRChem batched inference.

    This function enables efficient batched inference by submitting calculations
    concurrently via threads while the underlying Ray Serve batches requests
    to the GPU. This significantly improves throughput when applying the same
    calculation to many different structures.

    The InferenceBatcher is cached based on the model configuration, so repeated
    calls with the same model will reuse the same Ray Serve deployment.

    Parameters
    ----------
    func
        The function to map over atoms. Should accept `atoms` as first argument
        and `**calc_kwargs` with `method="fairchem"`. The function will receive
        `predict_unit` and `task_name` automatically injected into calc_kwargs.
    atoms_list
        List of ASE Atoms objects to process.
    fairchem_model
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
        Additional kwargs of the form key=list[...] that should be mapped over
        alongside atoms_list. Each list must have the same length as atoms_list.

    Returns
    -------
    list[Any]
        List of results from calling func for each atoms in atoms_list.

    Examples
    --------
    Basic usage with static_job:

    >>> from quacc.recipes.mlp.core import static_job
    >>> from quacc.wflow_tools.job_patterns import map_partitioned_lists_fairchembatch
    >>>
    >>> # Run batched static calculations
    >>> results = map_partitioned_lists_fairchembatch(
    ...     static_job,
    ...     atoms_list=my_atoms_list,
    ...     fairchem_model="uma-s-1",
    ...     task_name="omat",
    ... )

    With additional mapped kwargs:

    >>> results = map_partitioned_lists_fairchembatch(
    ...     my_custom_job,
    ...     atoms_list=my_atoms_list,
    ...     fairchem_model="uma-s-1",
    ...     task_name="omat",
    ...     custom_param=list_of_custom_values,  # mapped over
    ...     unmapped_kwargs={"shared_setting": value},  # same for all
    ... )
    """
    if not has_fairchem:
        raise ImportError(
            "fairchem must be installed to use map_partitioned_lists_fairchembatch. "
            "Run pip install fairchem-core."
        )

    from quacc.recipes.mlp._base import get_inference_batcher

    # Validate mapped_kwargs lengths
    if mapped_kwargs:
        all_lens = [len(v) for v in mapped_kwargs.values()]
        if not all(len(atoms_list) == le for le in all_lens):
            raise AssertionError(
                f"Inconsistent lengths: atoms_list has {len(atoms_list)} elements, "
                f"but mapped_kwargs have lengths {all_lens}"
            )

    # Get or create cached batcher
    batcher_kwargs = batcher_kwargs or {}
    batcher = get_inference_batcher(
        name_or_path=fairchem_model,
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
    for i, atoms in enumerate(atoms_list):
        # Build kwargs for this specific call
        call_kwargs = {**base_kwargs}
        for key, values in mapped_kwargs.items():
            call_kwargs[key] = values[i]

        future = batcher.executor.submit(raw_func, atoms, **call_kwargs)
        futures.append(future)

    # Wait for all results
    return [future.result() for future in futures]
