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


def make_volume_strains_from_bulk(
    atoms: Atoms,
    volume_strains: Sequence[float] = (
        -0.075,
        -0.06,
        -0.045,
        -0.03,
        -0.015,
        0.015,
        0.03,
        0.045,
        0.06,
        0.075,
    ),
) -> list[Atoms]:
    """
    Function to generate isotropically strained copies of a bulk atoms object.

    Parameters
    ----------
    atoms
        bulk atoms
    volume_strains
        Fractional changes in cell volume, e.g. -0.06 for a 6% compression. Each
        lattice vector is scaled by (1 + volume_strain) ** (1/3), and the
        fractional coordinates are kept fixed.

    Returns
    -------
    list[Atoms]
        The strained atoms objects, one per volume strain.
    """
    if any(strain <= -1 for strain in volume_strains):
        raise ValueError("Volume strains must be greater than -1.")

    strained_atoms = []
    for strain in volume_strains:
        new_atoms = atoms.copy()
        new_atoms.set_cell(atoms.get_cell() * (1 + strain) ** (1 / 3), scale_atoms=True)
        strained_atoms.append(new_atoms)
    return strained_atoms
