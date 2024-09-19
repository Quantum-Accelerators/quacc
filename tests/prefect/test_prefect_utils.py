from __future__ import annotations

import pytest
from prefect import flow, task

from quacc.wflow_tools.prefect_utils import (
    resolve_futures_to_results,
    resolve_futures_to_results_async,
)


def test_resolve_futures_to_results():
    @task
    def test_task():
        return {"test": 5}

    @flow
    def test_flow():
        future = test_task.submit()
        nested_future = {"nest": future}
        return resolve_futures_to_results(nested_future)

    result = test_flow()
    assert result["nest"]["test"] == 5


def test_resolve_futures_to_results_taskfail():
    @task
    def test_task():
        raise Exception
        return {"test": 5}

    @flow
    def test_flow():
        future = test_task.submit()
        nested_future = {"nest": future}
        return resolve_futures_to_results(nested_future)

    with pytest.raises(BaseException):  # noqa: B017
        test_flow()


async def test_resolve_futures_to_results_async():
    @task
    async def test_task():
        return {"test": 5}

    @flow
    async def test_flow():
        future = test_task.submit()
        nested_future = {"nest": future}
        return await resolve_futures_to_results_async(nested_future)

    result = await test_flow()
    assert result["nest"]["test"] == 5
