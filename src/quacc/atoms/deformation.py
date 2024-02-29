"""Utility functions for dealing with deformations."""

from __future__ import annotations

from collections.abc import Sequence

from ase import Atoms
from pymatgen.analysis.elasticity.strain import DeformedStructureSet
from pymatgen.io.ase import AseAtomsAdaptor


def make_deformations_from_bulk(
    atoms: Atoms,
    norm_strains: Sequence[float] = (-0.01, -0.005, 0.005, 0.01),
    shear_strains: Sequence[float] = (-0.06, -0.03, 0.03, 0.0),
    symmetry: bool = False,
) -> list[Atoms]:
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
    list[Atoms]
        All generated deformed structures
    """

    struct = AseAtomsAdaptor.get_structure(atoms)

    deformed_set = DeformedStructureSet(
        struct,
        norm_strains=norm_strains,
        shear_strains=shear_strains,
        symmetry=symmetry,
    )

    return [structure.to_ase_atoms() for structure in deformed_set]
