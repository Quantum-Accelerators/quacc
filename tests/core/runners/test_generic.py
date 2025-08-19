from __future__ import annotations

from subprocess import CalledProcessError

import pytest

from quacc import change_settings
from quacc.runners.generic import GenericRunner


def test_generic_runner(tmp_path, monkeypatch):
    """Test the `GenericRunner` class with a simple command."""
    gr = GenericRunner(command="echo Hello, World!")

    gr.filepaths["fd_out"] = "fd_out"
    gr.filepaths["fd_err"] = "fd_err"

    assert gr.tmpdir.exists()
    assert not gr.job_results_dir.exists()

    with change_settings({"GZIP_FILES": False}):
        gr.run_cmd()

    assert (gr.job_results_dir / "fd_out").read_text() == "Hello, World!\n"
    assert (gr.job_results_dir / "fd_err").read_text() == ""

    assert not gr.tmpdir.exists()
    assert gr.job_results_dir.exists()

    gr = GenericRunner(command="acommandthatdoesnotexistsurely")

    with pytest.raises(FileNotFoundError):
        gr.run_cmd()

    gr = GenericRunner(
        command='sh -c "echo VAR=${MY_VAR}; echo SHOULD_NOT=${SHOULD_NOT_EXIST:-missing}"',
        environment={"MY_VAR": "custom_value"},
    )

    gr.filepaths["fd_out"] = None
    gr.filepaths["fd_err"] = None

    results = gr.run_cmd()

    assert results.stdout == "VAR=custom_value\nSHOULD_NOT=missing\n"

    gr = GenericRunner(command="sh -c 'exit 1'")

    with pytest.raises(CalledProcessError) as exc_info:
        results = gr.run_cmd()

    e = exc_info.value
    assert e.returncode == 1
    assert "Command '['sh', '-c', 'exit 1']' returned non-zero exit status 1." in str(e)
