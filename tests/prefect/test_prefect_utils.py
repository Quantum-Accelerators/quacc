from __future__ import annotations

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
        return resolve_futures_to_results(future)

    result = test_flow()
    assert result["test"] == 5


async def test_resolve_futures_to_results_async():
    @task
    async def test_task():
        return {"test": 5}

    @flow
    async def test_flow():
        future = test_task.submit()
        return await resolve_futures_to_results_async(future)

    result = await test_flow()
    assert result["test"] == 5
