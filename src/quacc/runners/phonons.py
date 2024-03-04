from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc.runners.prep import calc_cleanup, calc_setup

try:
    import phonopy

    has_deps = phonopy is not None and find_spec("seekpath") is not None
except ImportError:
    has_deps = False

if TYPE_CHECKING:
    from numpy.typing import NDArray

    if phonopy:
        from phonopy import Phonopy


@requires(has_deps, "Phonopy or seekpath is not installed.")
def run_phonopy(
    phonon: Phonopy,
    forces: NDArray,
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
    atoms = None  # placeholder

    # Perform staging operations
    tmpdir, job_results_dir = calc_setup(atoms)

    # Run phonopy
    phonon.forces = forces
    phonon.produce_force_constants()
    phonon.run_mesh(with_eigenvectors=True)
    phonon.run_total_dos()
    phonon.run_thermal_properties(t_step=t_step, t_max=t_max, t_min=t_min)
    phonon.auto_band_structure(
        write_yaml=True, filename=Path(tmpdir, "phonopy_auto_band_structure.yaml")
    )
    phonon.save(Path(tmpdir, "phonopy.yaml"), settings={"force_constants": True})
    phonon.directory = job_results_dir

    # Perform cleanup operations
    calc_cleanup(atoms, tmpdir, job_results_dir)

    return phonon
