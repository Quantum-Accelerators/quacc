from __future__ import annotations

import pytest
from ase.build import bulk
from prefect import flow, task
from prefect.futures import resolve_futures_to_results

from quacc.recipes.emt.phonons import phonon_flow


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

    with pytest.raises(BaseException):  # noqa: B017, PT011
        test_flow()


async def test_resolve_futures_to_results_async():
    @task
    async def test_task():
        return {"test": 5}

    @flow
    async def test_flow():
        future = test_task.submit()
        nested_future = {"nest": future}
        return resolve_futures_to_results(nested_future)

    result = await test_flow()
    assert result["nest"]["test"] == 5


def test_phonon_flow():
    cu_bulk = bulk("Cu")
    phonon_flow(cu_bulk)
