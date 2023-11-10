"""Common workflows for phonons"""
from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import flow, job, subflow
from quacc.atoms.phonons import atoms_to_phonopy, phonopy_atoms_to_ase_atoms
from quacc.runners.ase import run_calc
from quacc.schemas.phonopy import summarize_phonopy

try:
    import phonopy
except ImportError:
    phonopy = None

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from numpy.typing import ArrayLike
    from phonopy import Phonopy

    from quacc.schemas._aliases.phonopy import PhononSchema


@flow
@requires(phonopy, "Phonopy must be installed. Run `pip install quacc[phonons]`")
def phonon_flow(
    atoms: Atoms,
    static_job: static_job,
    static_job_kwargs: dict[str, Any] | None = None,
    supercell_matrix: ArrayLike = ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
    atom_disp: float = 0.015,
    t_step: float = 10,
    t_min: float = 0,
    t_max: float = 1000,
    fields_to_store: dict[str, Any] = None,
) -> PhononSchema:
    """
    Calculate phonon properties.

    This module is adapted from `matcalc` (https://github.com/materialsvirtuallab/matcalc)

    Parameters
    ----------
    atoms
        Atoms object with calculator attached.
    static_job
        The static job to calculate the forces.
    static_job_kwargs
        The kwargs for `static_job`.
    supercell_matrix
        Supercell matrix to use. Defaults to 2x2x2 supercell.
    atom_disp
        Atomic displacement (A).
    t_step
        Temperature step (K).
    t_min
        Min temperature (K).
    t_max
        Max temperature (K).
    fields_to_store
        Fields to store in the database.

    Returns
    -------
    PhononSchema
        Dictionary of results from [quacc.schemas.phonopy.summarize_phonopy][]
    """

    static_job_kwargs = static_job_kwargs or {}

    @subflow
    def _phonopy_forces_subflow(atoms: Atoms) -> list:
        phonon = atoms_to_phonopy(atoms, supercell_matrix, atom_disp)
        supercells = [
            phonopy_atoms_to_ase_atoms(s) for s in phonon.supercells_with_displacements
        ]
        return [
            static_job(supercell, **static_job_kwargs)
            for supercell in supercells
            if supercell is not None
        ]

    @job
    def _phonopy_thermo_job(atoms: Atoms, forces: list["NDArray"]) -> PhononSchema:
        phonon = atoms_to_phonopy(atoms, supercell_matrix, atom_disp)
        phonon.forces = forces
        phonon.produce_force_constants()
        phonon.run_mesh()
        phonon.run_thermal_properties(t_step=t_step, t_max=t_max, t_min=t_min)

        return summarize_phonopy(
            phonon,
            calculator,
            input_atoms=input_atoms,
            additional_fields=fields_to_store,
        )

    forces = _phonopy_forces_subflow(atoms)
    return _phonopy_thermo_job(atoms, forces)
