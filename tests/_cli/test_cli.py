import pytest


@pytest.fixture
def test_yaml():
    from pathlib import Path

    return Path.cwd() / "test_quacc.yaml"


@pytest.fixture()
def runner():
    from typer.testing import CliRunner

    return CliRunner()


def setup_module():
    from quacc import SETTINGS

    SETTINGS.CONFIG_FILE = test_yaml


def teardown_module(default_settings):
    import os

    from quacc import SETTINGS

    if test_yaml.exists():
        os.remove(test_yaml)

    SETTINGS.CONFIG_FILE = default_settings.CONFIG_FILE


def test_version(runner):
    from quacc import __version__
    from quacc._cli.quacc import app

    response = runner.invoke(app, ["--version"])
    assert response.exit_code == 0
    assert __version__ in response.stdout


def test_help(runner):
    from quacc._cli.quacc import app

    response = runner.invoke(app, ["--help"])
    assert response.exit_code == 0


def test_set(runner):
    from quacc._cli.quacc import app

    response = runner.invoke(app, ["set", "WORKFLOW_ENGINE", "local"])
    assert response.exit_code == 0
    assert "local" in response.stdout
    val = None
    with open(test_yaml) as f:
        for line in f:
            if "WORKFLOW_ENGINE" in line:
                val = line.split(":")[-1].strip()
    assert val == "local"

    response = runner.invoke(app, ["set", "WORKFLOW_ENGINE", "covalent"])
    assert response.exit_code == 0
    assert "covalent" in response.stdout
    val = None
    with open(test_yaml) as f:
        for line in f:
            if "WORKFLOW_ENGINE" in line:
                val = line.split(":")[-1].strip()
    assert val == "covalent"

    response = runner.invoke(app, ["set", "VASP_PARALLEL_CMD", "dummy"])
    assert response.exit_code == 0
    assert "dummy" in response.stdout
    val = None
    with open(test_yaml) as f:
        for line in f:
            if "VASP_PARALLEL_CMD" in line:
                val = line.split(":")[-1].strip()
    assert val == "dummy"


def test_unset(runner):
    from quacc._cli.quacc import app

    response = runner.invoke(app, ["unset", "WORKFLOW_ENGINE"])
    assert response.exit_code == 0
    assert "WORKFLOW_ENGINE" in response.stdout
    lines = ""
    with open(test_yaml) as f:
        for _ in f:
            lines += ""
    assert "WORKFLOW_ENGINE" not in lines


def test_bad():
    from quacc._cli.quacc import app

    response = runner.invoke(app, ["set", "CONFIG_FILE", "here"])
    assert response.exit_code != 0
    response = runner.invoke(app, ["set", "bad", "dummy"])
    assert response.exit_code != 0
    response = runner.invoke(app, ["unset", "bad"])
    assert response.exit_code != 0
    response = runner.invoke(app, ["unset", "CONFIG_FILE"])
    assert response.exit_code != 0
