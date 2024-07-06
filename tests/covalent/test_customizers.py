from __future__ import annotations

import pytest

ct = pytest.importorskip("covalent")


from pathlib import Path

from covalent._workflow.lattice import Lattice

from quacc import flow, get_settings, job, strip_decorator, subflow
from quacc.wflow_tools.customizers import customize_funcs, redecorate


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

    assert hasattr(add, "electron_object")
    stripped_add = strip_decorator(add)
    assert stripped_add(1, 2) == 3
    assert not hasattr(stripped_add, "electron_object")

    assert isinstance(add2, Lattice)
    stripped_add2 = strip_decorator(add2)
    assert stripped_add2(1, 2) == 3
    assert not hasattr(stripped_add2, "electron_object")
    assert not isinstance(stripped_add2, Lattice)

    assert hasattr(add3, "electron_object")
    stripped_add3 = strip_decorator(add3)
    assert stripped_add3(1, 2) == 3
    assert not hasattr(stripped_add3, "electron_object")
    assert not isinstance(stripped_add3, Lattice)


def test_customize_jobs():
    @job
    def add(a, b):
        return a + b

    @job
    def mult(a, b):
        return a * b

    add_, mult_ = customize_funcs(
        ["add", "mult"], [add, mult], decorators={"add": job(executor="local")}
    )
    assert add_.electron_object.executor == "local"
    assert mult_.electron_object.executor is None

    add_, mult_ = customize_funcs(
        ["add", "mult"],
        [add, mult],
        decorators={"add": job(executor="local"), "mult": job(executor="local")},
    )
    assert add_.electron_object.executor == "local"
    assert mult_.electron_object.executor == "local"

    add_, mult_ = customize_funcs(
        ["add", "mult"], [add, mult], decorators={"all": job(executor="local")}
    )
    assert add_.electron_object.executor == "local"
    assert mult_.electron_object.executor == "local"

    add_, mult_ = customize_funcs(
        ["add", "mult"],
        [add, mult],
        decorators={"all": job(executor="local"), "mult": job(executor="dask")},
    )
    assert add_.electron_object.executor == "local"
    assert mult_.electron_object.executor == "dask"


def test_customize_flows():
    @job
    def add(a, b):
        return a + b

    @flow
    def workflow():
        return add(1, 2)

    workflow_ = customize_funcs(
        "workflow", workflow, decorators={"workflow": flow(executor="local")}
    )
    assert workflow_.metadata["executor"] == "local"


def test_customize_subflows():
    @job
    def add(a, b):
        return a + b

    @subflow
    def subworkflow():
        return add(1, 2)

    subworkflow_ = customize_funcs(
        "workflow", subworkflow, decorators={"workflow": subflow(executor="local")}
    )
    assert subworkflow_.electron_object.executor == "local"


def test_bad_customizers():
    @job
    def add(a, b):
        return a + b

    with pytest.raises(
        ValueError,
        match=r"Invalid decorator keys: \['bad'\]. Valid keys are: \['add'\]",
    ):
        customize_funcs("add", add, decorators={"bad": job(executor="local")})


def test_change_settings_redecorate_job(tmp_path_factory):
    tmp_dir1 = tmp_path_factory.mktemp("dir1")

    @job
    def write_file_job(name="job.txt"):
        with open(Path(get_settings().RESULTS_DIR, name), "w") as f:
            f.write("test file")

    write_file_job = redecorate(
        write_file_job, job(settings_swap={"RESULTS_DIR": tmp_dir1})
    )

    @flow
    def my_flow():
        return write_file_job()

    ct.get_result(ct.dispatch(my_flow)(), wait=True)
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

    # Test with redecorating a job in a flow
    ct.get_result(
        ct.dispatch(write_file_flow)(
            job_decorators={
                "write_file_job": job(settings_swap={"RESULTS_DIR": tmp_dir2})
            }
        ),
        wait=True,
    )
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

    @flow
    def my_flow():
        return write_file_job()

    ct.get_result(ct.dispatch(my_flow)(), wait=True)
    assert not Path(tmp_dir1 / "job.txt").exists()
    assert Path(tmp_dir2 / "job.txt").exists()
