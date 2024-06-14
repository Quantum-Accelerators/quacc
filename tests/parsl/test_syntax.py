from __future__ import annotations

import pytest

parsl = pytest.importorskip("parsl")


from pathlib import Path

from quacc import SETTINGS, flow, job, redecorate, strip_decorator, subflow
from quacc.wflow_tools.customizers import customize_funcs


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


def test_change_settings_redecorate(tmp_path_factory):
    tmp_dir1 = tmp_path_factory.mktemp("dir1")
    tmp_dir2 = tmp_path_factory.mktemp("dir2")

    @job
    def write_file_job(name="job.txt"):
        with open(Path(SETTINGS.RESULTS_DIR, name), "w") as f:
            f.write("test file")

    @flow
    def write_file_flow(name="flow.txt", job_decorators=None):
        write_file_job_ = customize_funcs(
            ["write_file_job"], [write_file_job], decorators=job_decorators
        )
        return write_file_job_(name=name)

    write_file_job = redecorate(
        write_file_job, job(settings_swap={"RESULTS_DIR": tmp_dir1})
    )

    write_file_job().result()
    assert Path(tmp_dir1 / "job.txt").exists()

    write_file_flow(
        job_decorators={"write_file_job": job(settings_swap={"RESULTS_DIR": tmp_dir2})}
    ).result()
    assert Path(tmp_dir2 / "flow.txt").exists()
