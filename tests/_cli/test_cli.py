import os
from pathlib import Path

from typer.testing import CliRunner

from quacc import SETTINGS, __version__
from quacc._cli.quacc import app

DEFAULT_SETTINGS = SETTINGS.copy()
test_yaml = Path.cwd() / "test_quacc.yaml"


def setup_module():
    SETTINGS.CONFIG_FILE = test_yaml


def teardown_module():
    SETTINGS = DEFAULT_SETTINGS
    if test_yaml.exists():
        os.remove(test_yaml)


runner = CliRunner()


def test_version():
    response = runner.invoke(app, ["--version"])
    assert response.exit_code == 0
    assert __version__ in response.stdout


def test_help():
    response = runner.invoke(app, ["--help"])
    assert response.exit_code == 0


def test_set():
    response = runner.invoke(app, ["set", "WORKFLOW_ENGINE", "local"])
    assert response.exit_code == 0
    assert "local" in response.stdout
    val = None
    with open(test_yaml, "r") as f:
        for line in f:
            if "WORKFLOW_ENGINE" in line:
                val = line.split(":")[-1].strip()
    assert val == "local"

    response = runner.invoke(app, ["set", "WORKFLOW_ENGINE", "covalent"])
    assert response.exit_code == 0
    assert "covalent" in response.stdout
    val = None
    with open(test_yaml, "r") as f:
        for line in f:
            if "WORKFLOW_ENGINE" in line:
                val = line.split(":")[-1].strip()
    assert val == "covalent"

    response = runner.invoke(app, ["set", "VASP_PARALLEL_CMD", "dummy"])
    assert response.exit_code == 0
    assert "dummy" in response.stdout
    val = None
    with open(test_yaml, "r") as f:
        for line in f:
            if "VASP_PARALLEL_CMD" in line:
                val = line.split(":")[-1].strip()
    assert val == "dummy"


def test_unset():
    response = runner.invoke(app, ["unset", "WORKFLOW_ENGINE"])
    assert response.exit_code == 0
    assert "WORKFLOW_ENGINE" in response.stdout
    lines = ""
    with open(test_yaml, "r") as f:
        for _ in f:
            lines += ""
    assert "WORKFLOW_ENGINE" not in lines


def test_bad():
    response = runner.invoke(app, ["set", "CONFIG_FILE", "here"])
    assert response.exit_code != 0
    response = runner.invoke(app, ["set", "bad", "dummy"])
    assert response.exit_code != 0
    response = runner.invoke(app, ["unset", "bad"])
    assert response.exit_code != 0
    response = runner.invoke(app, ["unset", "CONFIG_FILE"])
    assert response.exit_code != 0