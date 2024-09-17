from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any

from prefect.futures import PrefectFuture
from prefect.results import BaseResult
from prefect.utilities.annotations import quote
from prefect.utilities.collections import StopVisiting, visit_collection
from typing_extensions import TypeVar

if TYPE_CHECKING:
    from prefect.states import State

F = TypeVar("F")
R = TypeVar("R")


def resolve_futures_to_results(expr: PrefectFuture | Any) -> State | Any:
    """
    Given a Python built-in collection, recursively find `PrefectFutures` and build a
    new collection with the same structure with futures resolved to their final result.
    Resolving futures to their final result may wait for execution to complete.

    Unsupported object types will be returned without modification.

    This function is a trivial change from resolve_futures_to_states here:
    https://github.com/PrefectHQ/prefect/blob/main/src/prefect/futures.py
    """
    futures: set[PrefectFuture] = set()

    def _collect_futures(futures, expr, context):
        # Expressions inside quotes should not be traversed
        if isinstance(context.get("annotation"), quote):
            raise StopVisiting

        if isinstance(expr, PrefectFuture):
            futures.add(expr)

        return expr

    visit_collection(
        expr, visit_fn=partial(_collect_futures, futures), return_data=False, context={}
    )

    # if no futures were found, return the original expression
    if not futures:
        return expr

    # Get final states for each future
    results = []
    for future in futures:
        future.wait()
        result = future.state.result()
        if isinstance(result, BaseResult):
            result = result.get()
        results.append(result)

    states_by_future = dict(zip(futures, results, strict=False))

    def replace_futures_with_states(expr, context):
        # Expressions inside quotes should not be modified
        if isinstance(context.get("annotation"), quote):
            raise StopVisiting

        if isinstance(expr, PrefectFuture):
            return states_by_future[expr]
        else:
            return expr

    return visit_collection(
        expr, visit_fn=replace_futures_with_states, return_data=True, context={}
    )
