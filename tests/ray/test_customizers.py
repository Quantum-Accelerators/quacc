from __future__ import annotations

import pytest

ray = pytest.importorskip("ray")

from pathlib import Path

from quacc import flow, get_settings, job, redecorate, strip_decorator, subflow
from quacc.wflow_tools.customizers import customize_funcs, update_parameters


def test_strip_decorators():
    @job
    def add(a, b):
        return a + b

    @flow
    def add2(a, b):
        return a + b

    @subflow
    def add3(a, b):
        return a + b

    stripped_add = strip_decorator(add)
    assert stripped_add(1, 2) == 3

    stripped_add2 = strip_decorator(add2)
    assert stripped_add2(1, 2) == 3

    stripped_add3 = strip_decorator(add3)
    assert stripped_add3(1, 2) == 3


def test_customize_funcs(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b=1):
        return a + b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c, d=0):
        return [add(val, c) for val in vals]

    @flow
    def dynamic_workflow(a, b, c=1):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    @flow
    def test_dynamic_workflow(a, b, c=3):
        result1 = add(a, b)
        result2 = make_more(result1)
        return update_parameters(add_distributed, {"d": 1}, decorator="subflow")(
            result2, c
        )

    add_ = update_parameters(add, {"b": 3}, decorator="job")
    dynamic_workflow_ = update_parameters(dynamic_workflow, {"c": 4}, decorator="flow")
    assert (add_(1)).result() == 4
    # Note: d=1 update on inner subflow is not picked up because subflow runs on
    # driver and `add` already returns ObjectRefs; we just check it runs.
    assert (dynamic_workflow_(1, 2)).result() == [7, 7, 7]


def test_change_settings_redecorate_job(tmp_path_factory):
    tmp_dir1 = tmp_path_factory.mktemp("dir1")

    @job
    def write_file_job(name="job.txt"):
        with open(Path(get_settings().RESULTS_DIR, name), "w") as f:
            f.write("test file")

    write_file_job = redecorate(
        write_file_job, job(settings_swap={"RESULTS_DIR": tmp_dir1})
    )
    (write_file_job()).result()
    assert Path(tmp_dir1 / "job.txt").exists()


def test_change_settings_redecorate_flow(tmp_path_factory):
    tmp_dir2 = tmp_path_factory.mktemp("dir2")

    @job
    def write_file_job(name="job.txt"):
        with open(Path(get_settings().RESULTS_DIR, name), "w") as f:
            f.write("test file")

    @flow
    def write_file_flow(name="flow.txt", job_decorators=None):
        write_file_job_ = customize_funcs(
            ["write_file_job"], [write_file_job], decorators=job_decorators
        )
        return write_file_job_(name=name)

    (
        write_file_flow(
            job_decorators={
                "write_file_job": job(settings_swap={"RESULTS_DIR": tmp_dir2})
            }
        )
    ).result()
    assert Path(tmp_dir2 / "flow.txt").exists()


def test_double_change_settings_redecorate_job(tmp_path_factory):
    tmp_dir1 = tmp_path_factory.mktemp("dir1")
    tmp_dir2 = tmp_path_factory.mktemp("dir2")

    @job
    def write_file_job(name="job.txt"):
        with open(Path(get_settings().RESULTS_DIR, name), "w") as f:
            f.write("test file")

    write_file_job = redecorate(
        write_file_job, job(settings_swap={"RESULTS_DIR": tmp_dir1})
    )
    write_file_job = redecorate(
        write_file_job, job(settings_swap={"RESULTS_DIR": tmp_dir2})
    )
    (write_file_job()).result()
    assert not Path(tmp_dir1 / "job.txt").exists()
    assert Path(tmp_dir2 / "job.txt").exists()
