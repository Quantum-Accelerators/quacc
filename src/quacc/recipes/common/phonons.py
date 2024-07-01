"""Common workflows for phonons."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import job, subflow
from quacc.atoms.phonons import get_phonopy, phonopy_atoms_to_ase_atoms
from quacc.runners.phonons import PhonopyRunner
from quacc.schemas.phonons import summarize_phonopy

has_phonopy = bool(find_spec("phonopy"))
has_seekpath = bool(find_spec("seekpath"))

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc import Job
    from quacc.types import PhononSchema

    if has_phonopy:
        from phonopy import Phonopy


@subflow
@requires(has_phonopy, "Phonopy must be installed. Run `pip install quacc[phonons]`")
@requires(has_seekpath, "Seekpath must be installed. Run `pip install quacc[phonons]`")
def phonon_subflow(
    atoms: Atoms,
    force_job: Job,
    symprec: float = 1e-4,
    min_lengths: float | tuple[float, float, float] | None = 20.0,
    supercell_matrix: (
        tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]] | None
    ) = None,
    displacement: float = 0.01,
    t_step: float = 10,
    t_min: float = 0,
    t_max: float = 1000,
    phonopy_kwargs: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
) -> PhononSchema:
    """
    Calculate phonon properties.

    Parameters
    ----------
    atoms
        Atoms object with calculator attached.
    force_job
        The static job to calculate the forces.
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
        Additional fields to add to the output schema.

    Returns
    -------
    PhononSchema
        Dictionary of results from [quacc.schemas.phonons.summarize_phonopy][]
    """

    @subflow
    def _get_forces_subflow(supercells: list[Atoms]) -> list[dict]:
        return [
            force_job(supercell) for supercell in supercells if supercell is not None
        ]

    @job
    def _thermo_job(
        atoms: Atoms,
        phonopy: Phonopy,
        force_job_results: list[dict],
        t_step: float,
        t_min: float,
        t_max: float,
        additional_fields: dict[str, Any] | None,
    ) -> PhononSchema:
        parameters = force_job_results[-1].get("parameters")
        forces = [output["results"]["forces"] for output in force_job_results]
        phonopy_results = PhonopyRunner().run_phonopy(
            phonopy, forces, t_step=t_step, t_min=t_min, t_max=t_max
        )

        return summarize_phonopy(
            phonopy,
            atoms,
            phonopy_results.directory,
            parameters=parameters,
            additional_fields=additional_fields,
        )

    phonopy = get_phonopy(
        atoms,
        min_lengths=min_lengths,
        supercell_matrix=supercell_matrix,
        symprec=symprec,
        displacement=displacement,
        phonopy_kwargs=phonopy_kwargs,
    )
    supercells = [
        phonopy_atoms_to_ase_atoms(s) for s in phonopy.supercells_with_displacements
    ]
    force_job_results = _get_forces_subflow(supercells)
    return _thermo_job(
        atoms, phonopy, force_job_results, t_step, t_min, t_max, additional_fields
    )
