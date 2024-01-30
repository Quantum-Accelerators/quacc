"""Common workflows for phonons."""
from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import flow, job, subflow
from quacc.atoms.phonons import get_phonopy, phonopy_atoms_to_ase_atoms
from quacc.schemas.phonons import summarize_phonopy

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
    symprec: float = 1e-4,
    min_length: float | None = 15.0,
    displacement: float = 0.01,
    t_step: float = 10,
    t_min: float = 0,
    t_max: float = 1000,
    phonopy_kwargs: dict[str, Any] | None = None,
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
    min_length
        Minimum length of each lattice dimension (A).
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

    @subflow
    def _get_forces_subflow(atoms: Atoms) -> list[dict]:
        phonon = get_phonopy(
            atoms,
            min_length=min_length,
            symprec=symprec,
            displacement=displacement,
            phonopy_kwargs=phonopy_kwargs,
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
            atoms,
            min_length=min_length,
            symprec=symprec,
            displacement=displacement,
            phonopy_kwargs=phonopy_kwargs,
        )
        parameters = force_job_results[-1].get("parameters")
        forces = [output["results"]["forces"] for output in force_job_results]
        phonon.forces = forces
        phonon.produce_force_constants()
        phonon.run_mesh()
        phonon.run_thermal_properties(t_step=t_step, t_max=t_max, t_min=t_min)

        return summarize_phonopy(
            phonon, atoms, parameters=parameters, additional_fields=additional_fields
        )

    if relax_job is not None:
        atoms = relax_job(atoms)["atoms"]

    force_job_results = _get_forces_subflow(atoms)
    return _thermo_job(atoms, force_job_results)
