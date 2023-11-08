"""Phonon recipes for EMT"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase.calculators.emt import EMT

from quacc import flow
from quacc.recipes.common.phonons import common_phonon_flow

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
    **kwargs,
) -> PhononSchema:
    """
    Carry out a phonon calculation.

    Parameters
    ----------
    atoms
        Atoms object
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
    **kwargs
        Custom kwargs for the EMT calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `ase.calculators.emt.EMT` calculator.

        !!! Info "Calculator defaults"

            ```python
            {}
            ```

    Returns
    -------
    PhononSchema
        Dictionary of results from [quacc.schemas.phonopy.summarize_phonopy][]
    """

    calc = EMT(**kwargs)

    return common_phonon_flow(
        atoms,
        calc,
        supercell_matrix=supercell_matrix,
        atom_disp=atom_disp,
        t_step=t_step,
        t_min=t_min,
        t_max=t_max,
        fields_to_store={"name": "EMT Phonons"},
    )
