from __future__ import annotations

import os
import platform
from pathlib import Path

import pytest
from typer.testing import CliRunner

from quacc import __version__, _internally_set_settings
from quacc._cli.quacc import app

FILE_PATH = Path(__file__).parent
TEST_YAML = FILE_PATH / "test_quacc.yaml"


@pytest.fixture()
def runner():
    return CliRunner()


def setup_module():
    _internally_set_settings({"CONFIG_FILE": TEST_YAML})


def teardown_module():
    if TEST_YAML.exists():
        os.remove(TEST_YAML)
    _internally_set_settings(reset=True)


def test_version(runner):
    response = runner.invoke(app, ["--version"])
    assert response.exit_code == 0
    assert __version__ in response.stdout


def test_help(runner):
    response = runner.invoke(app, ["--help"])
    assert response.exit_code == 0


def test_set(runner):
    response = runner.invoke(app, ["set", "WORKFLOW_ENGINE", "None"])
    assert response.exit_code == 0
    assert "None" in response.stdout
    val = None
    with open(TEST_YAML) as f:
        for line in f:
            if "WORKFLOW_ENGINE" in line:
                val = line.split(":")[-1].strip()
    assert not val

    response = runner.invoke(app, ["set", "WORKFLOW_ENGINE", "covalent"])
    assert response.exit_code == 0
    assert "covalent" in response.stdout
    val = None
    with open(TEST_YAML) as f:
        for line in f:
            if "WORKFLOW_ENGINE" in line:
                val = line.split(":")[-1].strip()
    assert val == "covalent"

    response = runner.invoke(app, ["set", "VASP_PARALLEL_CMD", "dummy"])
    assert response.exit_code == 0
    assert "dummy" in response.stdout
    val = None
    with open(TEST_YAML) as f:
        for line in f:
            if "VASP_PARALLEL_CMD" in line:
                val = line.split(":")[-1].strip()
    assert val == "dummy"

    response = runner.invoke(app, ["set", "WORKFLOW_ENGINE", "invalid"])
    assert response.exit_code == 1

    response = runner.invoke(app, ["set", "GZIP_FILES", "True"])
    assert response.exit_code == 0
    assert "True" in response.stdout
    val = None
    with open(TEST_YAML) as f:
        for line in f:
            if "GZIP_FILES" in line:
                val = line.split(":")[-1].strip()
    assert val == "true"


def test_unset(runner):
    response = runner.invoke(app, ["unset", "WORKFLOW_ENGINE"])
    assert response.exit_code == 0
    assert "WORKFLOW_ENGINE" in response.stdout
    lines = ""
    with open(TEST_YAML) as f:
        for _ in f:
            lines += ""
    assert "WORKFLOW_ENGINE" not in lines


def test_info(runner):
    response = runner.invoke(app, ["info"])
    assert response.exit_code == 0
    assert __version__ in response.stdout
    assert platform.python_version() in response.stdout


def test_bad(runner):
    response = runner.invoke(app, ["set", "CONFIG_FILE", "here"])
    assert response.exit_code != 0
    response = runner.invoke(app, ["set", "bad", "dummy"])
    assert response.exit_code != 0
    response = runner.invoke(app, ["unset", "bad"])
    assert response.exit_code != 0
    response = runner.invoke(app, ["unset", "CONFIG_FILE"])
    assert response.exit_code != 0
