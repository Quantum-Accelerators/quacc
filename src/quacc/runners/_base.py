"""
Base class for runners.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc.runners.prep import calc_cleanup, calc_setup

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.utils.files import Filenames, SourceDirectory


class BaseRunner:
    """
    A base class for runners that sets up the calculation and cleans up the scratch directory.
    """

    def __init__(
        self,
        atoms: Atoms | None,
        copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    ) -> None:
        """
        Initialize the BaseRunner.

        Parameters
        ----------
        atoms
            The atoms object
        copy_files
            The files to copy to the scratch directory

        Returns
        -------
        None
        """
        self.atoms = atoms
        self.copy_files = copy_files
        self.tmpdir, self.job_results_dir = calc_setup(
            self.atoms, copy_files=self.copy_files
        )

    def cleanup(self) -> None:
        """
        Perform cleanup operations on the scratch directory.

        Returns
        -------
        None
        """
        calc_cleanup(self.atoms, self.tmpdir, self.job_results_dir)
