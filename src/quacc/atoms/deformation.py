"""Utility functions for dealing with deformations."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.analysis.elasticity.strain import DeformedStructureSet
from pymatgen.io.ase import AseAtomsAdaptor

if TYPE_CHECKING:
    from collections.abc import Sequence

    from ase.atoms import Atoms


def make_deformations_from_bulk(
    atoms: Atoms,
    norm_strains: Sequence[float] = (-0.01, -0.005, 0.005, 0.01),
    shear_strains: Sequence[float] = (-0.06, -0.03, 0.03, 0.06),
    symmetry: bool = False,
) -> DeformedStructureSet:
    """
    Function to generate deformed structures from a bulk atoms object.

    Parameters
    ----------
    atoms
        bulk atoms
    norm_strains
        strain values to apply to each normal mode.
    shear_strains
        strain values to apply to each shear mode.
    symmetry
        whether or not to use symmetry reduction

    Returns
    -------
    DeformedStructureSet
        A pymatgen DeformedStructureSet with information on the deformed
        structures and their strains, useful for fitting later.
    """
    struct = AseAtomsAdaptor.get_structure(atoms)  # type: ignore

    return DeformedStructureSet(
        struct,
        norm_strains=norm_strains,
        shear_strains=shear_strains,
        symmetry=symmetry,
    )
