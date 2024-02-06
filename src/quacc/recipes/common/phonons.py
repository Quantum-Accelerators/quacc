"""Common workflows for phonons."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from monty.dev import requires

from quacc import flow, job, subflow
from quacc.atoms.phonons import get_phonopy, phonopy_atoms_to_ase_atoms
from quacc.schemas.phonons import summarize_phonopy
from quacc.utils.dicts import recursive_dict_merge

try:
    import phonopy
except ImportError:
    phonopy = None

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc import Job
    from quacc.schemas._aliases.phonons import PhononSchema


@flow
@requires(phonopy, "Phonopy must be installed. Run `pip install quacc[phonons]`")
def phonon_flow(
    atoms: Atoms,
    force_job: Job,
    relax_job: Job | None = None,
    min_lengths: float | tuple[float, float, float] | None = 20.0,
    phonopy_kwargs: dict[str, Any] | None = None,
    force_constants_kwargs: dict[str, Any] | None = None,
    mesh_kwargs: dict[str, Any] | None = None,
    total_dos_kwargs: dict[str, Any] | None = None,
    thermal_properties_kwargs: dict[str, Any] | None = None,
    auto_band_structure_kwargs: dict[str, Any] | None = None,
    generate_displacements_kwargs: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
) -> PhononSchema:
    """
    Calculate phonon properties.

    This module is adapted from `matcalc` (https://github.com/materialsvirtuallab/matcalc)

    Parameters
    ----------
    atoms
        Atoms object with calculator attached.
    force_job
        The static job to calculate the forces.
    relax_job
        The job used to relax the structure before calculating the forces.
    symprec
        Precision for symmetry detection.
    min_lengths
        Minimum length of each lattice dimension (A).
    supercell_matrix
        The supercell matrix to use. If specified, it will override any
        value specified by `min_lengths`.
    displacement
        Atomic displacement (A).
    t_step
        Temperature step (K).
    t_min
        Min temperature (K).
    t_max
        Max temperature (K).
    phonopy_kwargs
        Additional kwargs to pass to the Phonopy class.
    additional_fields
        Additional fields to store in the database.

    Returns
    -------
    PhononSchema
        Dictionary of results from [quacc.schemas.phonons.summarize_phonopy][]
    """

    force_constants_kwargs = force_constants_kwargs or {}
    mesh_kwargs = mesh_kwargs or {}
    thermal_properties_kwargs = thermal_properties_kwargs or {}

    auto_band_structure_defaults = {
        "write_yaml": True,
        "filename": "phonopy_auto_band_structure.yaml",
    }

    auto_band_structure_kwargs = recursive_dict_merge(
        auto_band_structure_defaults, auto_band_structure_kwargs
    )

    if not set.intersection(
        set(thermal_properties_kwargs.keys()),
        {"t_step", "t_min", "t_max", "temperatures"},
    ):
        temperatures = np.linspace(0, 1000, 100, endpoint=True)
        np.append(temperatures, [273.15, 293.15, 298.15])
        thermal_properties_kwargs["temperatures"] = np.sort(temperatures)

    @subflow
    def _get_forces_subflow(atoms: Atoms) -> list[dict]:
        phonon = get_phonopy(
            atoms,
            min_lengths=min_lengths,
            phonopy_kwargs=phonopy_kwargs,
            generate_displacements_kwargs=generate_displacements_kwargs,
        )
        supercells = [
            phonopy_atoms_to_ase_atoms(s) for s in phonon.supercells_with_displacements
        ]
        return [
            force_job(supercell) for supercell in supercells if supercell is not None
        ]

    @job
    def _thermo_job(atoms: Atoms, force_job_results: list[dict]) -> PhononSchema:
        phonon = get_phonopy(
            atoms, min_lengths=min_lengths, phonopy_kwargs=phonopy_kwargs
        )
        parameters = force_job_results[-1].get("parameters")
        forces = [output["results"]["forces"] for output in force_job_results]
        phonon.forces = forces
        phonon.produce_force_constants(**force_constants_kwargs)
        phonon.run_mesh(**mesh_kwargs)
        phonon.run_total_dos(**total_dos_kwargs)
        phonon.run_thermal_properties(**thermal_properties_kwargs)
        phonon.auto_band_structure(
            **auto_band_structure_kwargs,
        )

        return summarize_phonopy(
            phonon, atoms, parameters=parameters, additional_fields=additional_fields
        )

    if relax_job is not None:
        atoms = relax_job(atoms)["atoms"]

    force_job_results = _get_forces_subflow(atoms)
    return _thermo_job(atoms, force_job_results)
