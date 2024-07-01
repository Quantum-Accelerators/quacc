"""
Base class for runners.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from quacc.runners.prep import calc_cleanup, calc_setup

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.types import Filenames, SourceDirectory


@dataclass
class BaseRunner:
    """
    A base class for runners that sets up the calculation and cleans up the scratch directory.

    Attributes
    ----------
    atoms
        Atoms object with calculator attached (or no Atoms object at all).
    copy_files
        Files to copy to runtime directory.
    """

    atoms: Atoms | None = None
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None

    def setup(self) -> None:
        """
        Perform setup operations on the runtime directory.

        Returns
        -------
        None
        """
        self.tmpdir, self.job_results_dir = calc_setup(
            self.atoms, copy_files=self.copy_files
        )

    def cleanup(self) -> None:
        """
        Perform cleanup operations on the runtime directory.

        Returns
        -------
        None
        """
        calc_cleanup(self.atoms, self.tmpdir, self.job_results_dir)
