from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase.units import _c as c
from monty.dev import requires

from quacc.runners.prep import calc_cleanup, calc_setup
from quacc.utils.dicts import recursive_dict_merge

try:
    import phonopy

    has_deps = phonopy is not None and find_spec("seekpath") is not None
except ImportError:
    has_deps = False

if TYPE_CHECKING:
    from typing import Any

    from numpy.typing import NDArray

    if phonopy:
        from phonopy import Phonopy


@requires(has_deps, "Phonopy or seekpath is not installed.")
def run_phonopy(
    phonon: Phonopy,
    forces: NDArray,
    is_fixed: bool = False,
    force_constants_kwargs: dict[str, Any] | None = None,
    mesh_kwargs: dict[str, Any] | None = None,
    total_dos_kwargs: dict[str, Any] | None = None,
    thermal_properties_kwargs: dict[str, Any] | None = None,
    auto_band_structure_kwargs: dict[str, Any] | None = None,
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
    force_constants_kwargs
        Additional kwargs to pass to the `produce_force_constants` method.
    mesh_kwargs
        Additional kwargs to pass to the `run_mesh` method.
    total_dos_kwargs
        Additional kwargs to pass to the `run_total_dos` method.
    thermal_properties_kwargs
        Additional kwargs to pass to the `run_thermal_properties` method.
    auto_band_structure_kwargs
        Additional kwargs to pass to the `auto_band_structure` method.

    Returns
    -------
    Phonopy
        The phonopy object with the results.
    """
    atoms = None  # placeholder

    # Perform staging operations
    tmpdir, job_results_dir = calc_setup(atoms)

    force_constants_kwargs = force_constants_kwargs or {}
    thermal_properties_kwargs = thermal_properties_kwargs or {}

    auto_band_structure_defaults = {
        "write_yaml": True,
        "filename": Path(tmpdir, "phonopy_auto_band_structure.yaml"),
    }

    auto_band_structure_kwargs = recursive_dict_merge(
        auto_band_structure_defaults, auto_band_structure_kwargs
    )

    total_dos_defaults = {"freq_pitch": c * 100 * 1.0e-12}
    total_dos_kwargs = recursive_dict_merge(total_dos_defaults, total_dos_kwargs)

    if not set.intersection(
        set(thermal_properties_kwargs.keys()),
        {"t_step", "t_min", "t_max", "temperatures"},
    ):
        temperatures = np.append(
            np.linspace(0, 1000, 101, endpoint=True), [273.15, 293.15, 298.15]
        )
        thermal_properties_kwargs["temperatures"] = np.sort(temperatures)

    mesh_defaults = {"with_eigenvectors": True}
    mesh_kwargs = recursive_dict_merge(mesh_defaults, mesh_kwargs)

    # Run phonopy
    phonon.forces = forces
    phonon.produce_force_constants(**force_constants_kwargs)

    if is_fixed:
        phonon.symmetrize_force_constants()
        phonon.symmetrize_force_constants_by_space_group()

    phonon.run_mesh(**mesh_kwargs)
    phonon.run_total_dos(**total_dos_kwargs)
    phonon.run_thermal_properties(**thermal_properties_kwargs)
    phonon.auto_band_structure(**auto_band_structure_kwargs)
    phonon.save(Path(tmpdir, "phonopy.yaml"), settings={"force_constants": True})
    phonon.directory = job_results_dir

    # Perform cleanup operations
    calc_cleanup(atoms, tmpdir, job_results_dir)

    return phonon
