"""Utility functions for thermochemistry."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ase import units
from ase.thermochemistry import IdealGasThermo

from quacc.schemas.atoms import atoms_to_metadata

if TYPE_CHECKING:
    from ase.atoms import Atoms


def run_ideal_gas(
    atoms: Atoms,
    vib_freqs: list[float | complex],
    energy: float = 0.0,
    spin_multiplicity: int | None = None,
) -> IdealGasThermo:
    """
    Create an IdealGasThermo object for a molecule from a given vibrational analysis.
    This is for free gases only and will not be valid for solids or adsorbates on
    surfaces. Any imaginary vibrational modes after the 3N-5/3N-6 cut will simply be
    ignored.

    Parameters
    ----------
    atoms
        The Atoms object associated with the vibrational analysis.
    vib_freqs
        The list of vibrations to use in cm^-1, typically obtained from
        Vibrations.get_frequencies().
    energy
        Potential energy in eV. If 0 eV, then the thermochemical correction is
        computed.
    spin_multiplicity
        The spin multiplicity (2S+1). If None, this will be determined
        automatically from the attached magnetic moments.

    Returns
    -------
    IdealGasThermo object
    """

    # Switch off PBC since this is only for molecules
    atoms.set_pbc(False)

    # Ensure all negative modes are made complex
    for i, f in enumerate(vib_freqs):
        if not isinstance(f, complex) and f < 0:
            vib_freqs[i] = complex(0 - f * 1j)

    # Convert vibrational frequencies to energies
    vib_energies = [f * units.invcm for f in vib_freqs]

    # Get the spin from the Atoms object.
    if spin_multiplicity:
        spin = (spin_multiplicity - 1) / 2
    elif (
        getattr(atoms, "calc", None) is not None
        and getattr(atoms.calc, "results", None) is not None
        and atoms.calc.results.get("magmom", None) is not None
    ):
        spin = round(atoms.calc.results["magmom"]) / 2
    elif (
        getattr(atoms, "calc", None) is not None
        and getattr(atoms.calc, "results", None) is not None
        and atoms.calc.results.get("magmoms", None) is not None
    ):
        spin = round(np.sum(atoms.calc.results["magmoms"])) / 2
    elif atoms.has("initial_magmoms"):
        spin = round(np.sum(atoms.get_initial_magnetic_moments())) / 2
    else:
        spin = 0

    # Get symmetry for later use
    natoms = len(atoms)
    metadata = atoms_to_metadata(atoms)

    # Get the geometry
    if natoms == 1:
        geometry = "monatomic"
    elif metadata["symmetry"]["linear"]:
        geometry = "linear"
    else:
        geometry = "nonlinear"

    return IdealGasThermo(
        vib_energies,
        geometry,
        potentialenergy=energy,
        atoms=atoms,
        symmetrynumber=metadata["symmetry"]["rotation_number"],
        spin=spin,
        ignore_imag_modes=True,
    )
