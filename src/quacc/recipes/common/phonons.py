"""Common workflows for phonons"""
from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure

from quacc import flow, job, subflow
from quacc.recipes.common.core import force_job
from quacc.schemas.phonopy import summarize_phonopy

try:
    import phonopy
except ImportError:
    phonopy = None

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from ase.calculators.calculator import Calculator
    from numpy.typing import ArrayLike, NDArray
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    from quacc.schemas._aliases.phonopy import PhononSchema


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

    @subflow
    def _calc_forces_distributed(supercells: list[Atoms]) -> list[NDArray]:
        return [
            force_job(supercell, calculator)
            for supercell in supercells
            if supercell is not None
        ]

    @job
    def _run_phonons(atoms: Atoms) -> PhononSchema:
        phonon = _atoms_to_phonopy(atoms, supercell_matrix, atom_disp)
        supercells = [
            _phonopy_atoms_to_ase_atoms(s) for s in phonon.supercells_with_displacements
        ]
        forces = _calc_forces_distributed(supercells)
        phonon.forces = forces
        phonon.produce_force_constants()
        phonon.run_mesh()
        phonon.run_thermal_properties(t_step=t_step, t_max=t_max, t_min=t_min)
        return summarize_phonopy(
            phonon,
            calculator,
            input_atoms=atoms,
            additional_fields=fields_to_store,
        )

    return _run_phonons(atoms)


def _atoms_to_phonopy(
    atoms: Atoms, supercell_matrix: ArrayLike, atom_disp: float
) -> Phonopy:
    """
    Convert an ASE atoms object to a phonopy object with displacements
    generated.

    Parameters
    ----------
    atoms
        ASE atoms object
    supercell_matrix
        Supercell matrix to use.
    atom_disp
        Atomic displacement (A).

    Returns
    -------
    Phonopy
        Phonopy object
    """
    structure = AseAtomsAdaptor().get_structure(atoms)
    phonopy_atoms = get_phonopy_structure(structure)
    phonon = phonopy.Phonopy(phonopy_atoms, supercell_matrix)
    phonon.generate_displacements(distance=atom_disp)
    return phonon


def _phonopy_atoms_to_ase_atoms(phonpy_atoms: PhonopyAtoms) -> Atoms:
    """
    Convert a phonopy atoms object to an ASE atoms object.

    Parameters
    ----------
    phonpy_atoms
        Phonopy atoms object

    Returns
    -------
    Atoms
        ASE atoms object
    """
    pmg_structure = get_pmg_structure(phonpy_atoms)
    return AseAtomsAdaptor().get_atoms(pmg_structure)
