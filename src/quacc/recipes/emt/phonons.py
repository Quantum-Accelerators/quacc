"""Phonon recipes for EMT"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT

from quacc import flow
from quacc.recipes.common.phonons import run_phonons
from quacc.schemas.phonopy import summarize_phonopy

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms
    from numpy.typing import ArrayLike

    from quacc.recipes.phonons import PhononSchema


@flow
def phonon_flow(
    atoms: Atoms,
    supercell_matrix: ArrayLike = ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
    atom_disp: float = 0.015,
    t_step: float = 10,
    t_min: float = 0,
    t_max: float = 1000,
    calc_swaps: dict[str, Any] | None = None,
) -> PhononSchema:
    """
    Carry out a phonon calculation.

    Parameters
    ----------
    atoms
        Atoms object
    calc_swaps
        Dictionary of custom kwargs for the EMT calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.emt.EMT` calculator.

        !!! Info "Calculator defaults"

            ```python
            {}
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    PhononSchema
        Dictionary of results from [quacc.schemas.phonopy.summarize_phonopy][]
    """

    calc_swaps = calc_swaps or {}
    atoms.calc = EMT(**calc_swaps)

    phonon = run_phonons(
        atoms,
        supercell_matrix=supercell_matrix,
        atom_disp=atom_disp,
        t_step=t_step,
        t_min=t_min,
        t_max=t_max,
    )
    return summarize_phonopy(
        phonon,
        input_atoms=atoms,
        additional_fields={"name": "EMT Phonons"},
    )
