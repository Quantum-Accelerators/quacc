from __future__ import annotations

import pickle
from functools import partial

import pytest

ray = pytest.importorskip("ray")

from quacc import flow, job, subflow
from quacc.wflow_tools.decorators import (
    RayFuture,
    _unwrap_ray_future,
    _wrap_partial_for_ray,
)


def test_unwrap_ray_future_passthrough():
    assert _unwrap_ray_future(5) == 5
    assert _unwrap_ray_future("abc") == "abc"
    assert _unwrap_ray_future(None) is None


def test_unwrap_ray_future_containers():
    @job
    def add(a, b):
        return a + b

    fut1 = add(1, 2)
    fut2 = add(3, 4)

    assert isinstance(fut1, RayFuture)

    out_list = _unwrap_ray_future([fut1, 5, fut2])
    assert out_list[0] is fut1._ref
    assert out_list[1] == 5
    assert out_list[2] is fut2._ref

    out_tuple = _unwrap_ray_future((fut1, "x"))
    assert out_tuple[0] is fut1._ref
    assert out_tuple[1] == "x"

    out_dict = _unwrap_ray_future({"a": fut1, "b": 7})
    assert out_dict["a"] is fut1._ref
    assert out_dict["b"] == 7

    nested = _unwrap_ray_future({"k": [fut1, (fut2, 1)]})
    assert nested["k"][0] is fut1._ref
    assert nested["k"][1][0] is fut2._ref


def test_unwrap_ray_future_subclasses_passthrough():
    class MyList(list):
        pass

    class MyDict(dict):
        pass

    ml = MyList([1, 2])
    md = MyDict({"a": 1})
    assert _unwrap_ray_future(ml) is ml
    assert _unwrap_ray_future(md) is md


def test_wrap_partial_for_ray_noop_for_plain():
    def f(a, b):
        return a + b

    assert _wrap_partial_for_ray(f) is f


def test_wrap_partial_for_ray_with_partial():
    def f(a, b):
        return a + b

    p = partial(f, 10)
    wrapped = _wrap_partial_for_ray(p)
    assert wrapped is not p
    assert callable(wrapped)
    assert wrapped(5) == 15
    assert wrapped.__name__ == f.__name__


def test_rayfuture_reduce_pickles():
    @job
    def make_value():
        return 42

    fut = make_value()
    assert isinstance(fut, RayFuture)
    data = pickle.dumps(fut)
    restored = pickle.loads(data)
    assert isinstance(restored, RayFuture)
    assert ray.get(restored._ref) == 42


def test_rayfuture_getitem(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def make_dict():
        return {"a": 1, "b": [10, 20]}

    fut = make_dict()
    a_fut = fut["a"]
    assert isinstance(a_fut, RayFuture)
    assert a_fut.result() == 1

    b_fut = fut["b"]
    assert b_fut.result() == [10, 20]
    assert b_fut[1].result() == 20


def test_subflow_returns_tuple(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @subflow
    def make_tuple(x):
        return (add(x, 1), add(x, 2))

    @flow
    def wf(x):
        return make_tuple(x)

    result = wf(5).result()
    assert result == (6, 7)


def test_subflow_returns_dict(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @subflow
    def make_dict(x):
        return {"first": add(x, 1), "second": add(x, 2)}

    @flow
    def wf(x):
        return make_dict(x)

    result = wf(5).result()
    assert result == {"first": 6, "second": 7}


def test_subflow_returns_scalar(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @subflow
    def single(x):
        return add(x, 100)

    @flow
    def wf(x):
        return single(x)

    assert wf(5).result() == 105


def test_subflow_with_kwargs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @subflow(num_cpus=1)
    def make_pair(x):
        return [add(x, 1), add(x, 2)]

    @flow
    def wf(x):
        return make_pair(x)

    assert wf(3).result() == [4, 5]


def test_job_with_kwargs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job(num_cpus=1)
    def add(a, b):
        return a + b

    assert add(2, 3).result() == 5
