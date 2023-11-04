"""Common workflows for phonons"""
from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure

from quacc import flow, subflow
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
    from numpy.typing import ArrayLike
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
        Atomic displacement.
    t_step
        Temperature step.
    t_min
        Min temperature.
    t_max
        Max temperature.
    fields_to_store
        Fields to store in the database.

    Returns
    -------
    PhononSchema
        Dictionary of results from [quacc.schemas.phonopy.summarize_phonopy][]
    """

    @subflow
    def _calc_phonons_distributed(atoms: Atoms) -> PhononSchema:
        structure = AseAtomsAdaptor().get_structure(atoms)

        phonopy_atoms = get_phonopy_structure(structure)
        phonon = phonopy.Phonopy(phonopy_atoms, supercell_matrix)
        phonon.generate_displacements(distance=atom_disp)
        supercells = phonon.supercells_with_displacements

        phonon.forces = [
            force_job(_get_atoms_from_phonopy(supercell), calculator)
            for supercell in supercells
            if supercell is not None
        ]
        phonon.produce_force_constants()
        phonon.run_mesh()
        phonon.run_thermal_properties(t_step=t_step, t_max=t_max, t_min=t_min)
        return summarize_phonopy(
            phonon,
            calculator,
            input_atoms=atoms,
            additional_fields=fields_to_store,
        )

    return _calc_phonons_distributed(atoms)


def _get_atoms_from_phonopy(phonpy_atoms: PhonopyAtoms) -> Atoms:
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
