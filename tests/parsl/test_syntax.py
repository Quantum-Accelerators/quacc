from __future__ import annotations

import pytest

parsl = pytest.importorskip("parsl")

import os
from pathlib import Path

from quacc import SETTINGS, change_settings_wf, flow, job, strip_decorator, subflow


def test_parsl_decorators(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    @job
    def make_more(val):
        return [val] * 3

    @subflow
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @subflow
    def add_distributed2(vals, c, op):
        return [op(val, c) for val in vals]

    @flow
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    @flow
    def dynamic_workflow2(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed2(result2, c, add)

    @flow
    def dynamic_workflow3(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        result3 = add_distributed(result2, c)
        result4 = add_distributed(result3, c)
        return add(result4[0], c)

    assert add(1, 2).result() == 3
    assert mult(1, 2).result() == 2
    assert workflow(1, 2, 3).result() == 9
    assert dynamic_workflow(1, 2, 3).result() == [6, 6, 6]
    assert dynamic_workflow2(1, 2, 3).result() == [6, 6, 6]
    assert dynamic_workflow3(1, 2, 3).result() == 12


def test_parsl_decorators_args(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    @job()
    def add(a, b):
        return a + b

    @job()
    def mult(a, b):
        return a * b

    @job()
    def make_more(val):
        return [val] * 3

    @subflow()
    def add_distributed(vals, c):
        return [add(val, c) for val in vals]

    @flow()
    def workflow(a, b, c):
        return mult(add(a, b), c)

    @flow()
    def dynamic_workflow(a, b, c):
        result1 = add(a, b)
        result2 = make_more(result1)
        return add_distributed(result2, c)

    assert add(1, 2).result() == 3
    assert mult(1, 2).result() == 2
    assert workflow(1, 2, 3).result() == 9
    assert dynamic_workflow(1, 2, 3).result() == [6, 6, 6]


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


def test_special_params(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)

    @job(walltime=30, parsl_resource_specification={})
    def add(a, b):
        return a + b

    @subflow(walltime=20, parsl_resource_specification={})
    def add2(a, b):
        return [add(i, i + 2) for i in range(a, b + 5)]

    assert add(1, 2).result() == 3
    assert add2(1, 2).result() == [4, 6, 8, 10, 12, 14]


def test_change_settings_wf(tmp_path_factory):
    @job
    def write_file(name="job"):
        with open(Path(f"{SETTINGS.RESULTS_DIR}/{name}.txt"), "w") as f:
            f.write("test file")

    @flow
    def write_file2():
        with open(Path(f"{SETTINGS.RESULTS_DIR}/flow.txt"), "w") as f:
            f.write("test file")

    @subflow
    def write_file3():
        return [write_file(name="subflow")]

    tmp_dir = tmp_path_factory.mktemp("dir1")
    tmp_dir2 = tmp_path_factory.mktemp("dir2")
    tmp_dir3 = tmp_path_factory.mktemp("dir3")

    write_file_new = change_settings_wf(write_file, {"RESULTS_DIR": tmp_dir}, job)
    write_file2_new = change_settings_wf(write_file2, {"RESULTS_DIR": tmp_dir2}, flow)
    write_file3_new = change_settings_wf(
        write_file3, {"RESULTS_DIR": tmp_dir3}, subflow
    )

    write_file_new().result()
    write_file2_new()
    write_file3_new().result()

    assert os.path.exists(tmp_dir / "job.txt")
    assert os.path.exists(tmp_dir2 / "flow.txt")
    assert os.path.exists(tmp_dir3 / "subflow.txt")
