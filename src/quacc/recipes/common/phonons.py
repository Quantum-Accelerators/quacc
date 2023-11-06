"""Common workflows for phonons"""
from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import flow, job, subflow
from quacc.atoms.phonons import atoms_to_phonopy, phonopy_atoms_to_ase_atoms
from quacc.recipes.common.core import force_job
from quacc.schemas.phonopy import summarize_phonopy

try:
    import phonopy
except ImportError:
    phonopy = None

if TYPE_CHECKING:
    from typing import Any, TypedDict

    from ase.atoms import Atoms
    from ase.calculators.calculator import Calculator
    from numpy.typing import ArrayLike, NDArray
    from phonopy import Phonopy

    from quacc.schemas._aliases.phonopy import PhononSchema

    class ForceReturn(TypedDict):
        phonon: Phonopy
        forces: list[NDArray]


@flow
@requires(phonopy, "Phonopy must be installed. Run `pip install quacc[phonons]`")
def phonon_flow(
    atoms: Atoms,
    calculator: Calculator,
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
    calculator
        Calculator to use.
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
    fields_to_store = fields_to_store or {}

    @subflow
    def _force_job_distributed(atoms: Atoms) -> ForceReturn:
        phonon = atoms_to_phonopy(atoms, supercell_matrix, atom_disp)
        supercells = [
            phonopy_atoms_to_ase_atoms(s) for s in phonon.supercells_with_displacements
        ]
        return {
            "phonon": phonon,
            "forces": [
                force_job(supercell, calculator)
                for supercell in supercells
                if supercell is not None
            ],
        }

    @job
    def _thermo_job(force_return: ForceReturn, input_atoms: Atoms) -> PhononSchema:
        phonon = force_return["phonon"]
        phonon.forces = force_return["forces"]
        phonon.produce_force_constants()
        phonon.run_mesh()
        phonon.run_thermal_properties(t_step=t_step, t_max=t_max, t_min=t_min)

        return summarize_phonopy(
            phonon,
            calculator,
            input_atoms=input_atoms,
            additional_fields=fields_to_store,
        )

    force_return = _force_job_distributed(atoms)

    return _thermo_job(force_return, atoms)
