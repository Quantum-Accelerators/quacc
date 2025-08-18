"""Class to run generic commands in a subprocess."""

from __future__ import annotations

import os
from contextlib import ExitStack
from pathlib import Path
from shlex import split
from subprocess import PIPE, CompletedProcess, run
from typing import TYPE_CHECKING, ClassVar, Final

from quacc.runners._base import BaseRunner

if TYPE_CHECKING:
    from quacc.types import Filenames, SourceDirectory


class GenericRunner(BaseRunner):
    """
    A class to run generic (IO) commands in a subprocess. Inherits from BaseRunner, which handles setup and cleanup of the calculation.
    """

    filepaths: ClassVar[dict[str, SourceDirectory | None]] = {
        "fd_out": None,
        "fd_err": None,
    }

    def __init__(
        self,
        command: str,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
        environment: dict[str, str] | None = None,
    ) -> None:
        """
        Initialize the `GenericRunner` with the command and optional copy files and environment variables.
        
        Parameters
        ----------
        command
            The command to run in the subprocess.
        copy_files
            Files to copy to the runtime directory.
        environment
            Environment variables to set for the subprocess. If None, the current environment is used.

        Returns
        -------
        None
        """
        self.command: Final[list[str]] = split(command)
        self.environment = environment

        super().__init__(copy_files=copy_files)

        self.setup()

    def run_cmd(self) -> CompletedProcess:
        """
        Run a command in a subprocess.

        Returns
        -------
        CompletedProcess
            The result of the subprocess execution.
        """
        with ExitStack() as stack:
            files = {
                name: stack.enter_context(Path(self.tmpdir, path).open("w"))
                for name, path in self.filepaths.items()
                if path is not None
            }

            cmd_results = run(
                self.command,
                cwd=self.tmpdir,
                shell=False,
                check=True,
                env=self.environment if self.environment is not None else os.environ,
                stdout=files.get("fd_out", PIPE),
                stderr=files.get("fd_err", PIPE),
                text=True,
            )

        self.cleanup()
        return cmd_results
