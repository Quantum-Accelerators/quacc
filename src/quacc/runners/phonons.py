from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc.runners.ase import BaseRunner

has_phonopy = bool(find_spec("phonopy"))
has_seekpath = bool(find_spec("seekpath"))


if TYPE_CHECKING:
    from numpy.typing import NDArray

    if has_phonopy:
        from phonopy import Phonopy


class PhonopyRunner(BaseRunner):
    def __init__(self) -> None:
        """
        Initialize the PhonopyRunner.

        Returns
        -------
        None
        """
        self.setup()

    @requires(has_phonopy, "Phonopy is not installed.")
    @requires(has_seekpath, "Seekpath is not installed")
    def run_phonopy(
        self,
        phonon: Phonopy,
        forces: NDArray,
        symmetrize: bool = False,
        t_step: float = 10,
        t_min: float = 0,
        t_max: float = 1000,
    ) -> Phonopy:
        """
        Run a phonopy calculation in a temporary directory and
        copy the results to the job results directory.

        Parameters
        ----------
        phonon
            Phonopy object
        forces
            Forces on the atoms
        symmetrize
            Whether to symmetrize the force constants
        t_step
            Temperature step
        t_min
            Minimum temperature
        t_max
            Maximum temperature

        Returns
        -------
        Phonopy
            The phonopy object with the results.
        """

        # Run phonopy
        phonon.forces = forces
        phonon.produce_force_constants()

        if symmetrize:
            phonon.symmetrize_force_constants()
            phonon.symmetrize_force_constants_by_space_group()

        phonon.run_mesh(with_eigenvectors=True)
        phonon.run_total_dos()
        phonon.run_thermal_properties(t_step=t_step, t_max=t_max, t_min=t_min)
        phonon.auto_band_structure(
            write_yaml=True,
            filename=Path(self.tmpdir, "phonopy_auto_band_structure.yaml"),
        )
        phonon.save(
            Path(self.tmpdir, "phonopy.yaml"), settings={"force_constants": True}
        )
        phonon.directory = self.job_results_dir

        # Perform cleanup operations
        self.cleanup()

        return phonon
