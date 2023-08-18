from click.testing import CliRunner

from quacc._cli.cli import cli


def test_help():
    response = CliRunner().invoke(cli).output
    assert "Show this message and exit" in response
