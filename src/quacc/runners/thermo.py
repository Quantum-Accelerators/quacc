"""Utility functions for thermochemistry."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from ase import units
from ase.thermochemistry import IdealGasThermo, HarmonicThermo
from emmet.core.symmetry import PointGroupData
from pymatgen.io.ase import AseAtomsAdaptor

if TYPE_CHECKING:
    from ase.atoms import Atoms


@dataclass
class ThermoRunner:
    """
    A runner for thermochemistry calculations.

    Attributes
    ----------
    atoms
        The Atoms object associated with the vibrational analysis.
    vib_freqs
        The list of vibrations to use in cm^-1, typically obtained from
        Vibrations.get_frequencies().
    energy
        Potential energy in eV. If 0 eV, then the thermochemical correction is
        computed.
    """

    atoms: Atoms
    vib_freqs: list[float | complex]
    energy: float = 0.0

    def run_ideal_gas(self, spin_multiplicity: int | None = None) -> IdealGasThermo:
        """
        Create an IdealGasThermo object for a molecule from a given vibrational analysis.
        This is for free gases only and will not be valid for solids or adsorbates on
        surfaces. Any imaginary vibrational modes after the 3N-5/3N-6 cut will simply be
        ignored.

        Parameters
        ----------
        spin_multiplicity
            The spin multiplicity (2S+1). If None, this will be determined
            automatically from the attached magnetic moments.

        Returns
        -------
        IdealGasThermo object
        """
        # Ensure all negative modes are made complex
        for i, f in enumerate(self.vib_freqs):
            if not isinstance(f, complex) and f < 0:
                self.vib_freqs[i] = complex(0 - f * 1j)

        # Convert vibrational frequencies to energies
        vib_energies = [f * units.invcm for f in self.vib_freqs]

        # Get the spin from the Atoms object.
        if spin_multiplicity:
            spin = (spin_multiplicity - 1) / 2
        elif (
            getattr(self.atoms, "calc", None) is not None
            and getattr(self.atoms.calc, "results", None) is not None
            and self.atoms.calc.results.get("magmom", None) is not None
        ):
            spin = round(self.atoms.calc.results["magmom"]) / 2
        elif (
            getattr(self.atoms, "calc", None) is not None
            and getattr(self.atoms.calc, "results", None) is not None
            and self.atoms.calc.results.get("magmoms", None) is not None
        ):
            spin = round(np.sum(self.atoms.calc.results["magmoms"])) / 2
        elif self.atoms.has("initial_magmoms"):
            spin = round(np.sum(self.atoms.get_initial_magnetic_moments())) / 2
        else:
            spin = 0

        # Get symmetry for later use
        natoms = len(self.atoms)
        mol = AseAtomsAdaptor().get_molecule(self.atoms, charge_spin_check=False)
        point_group_data = PointGroupData().from_molecule(mol)

        # Get the geometry
        if natoms == 1:
            geometry = "monatomic"
        elif point_group_data.linear:
            geometry = "linear"
        else:
            geometry = "nonlinear"

        return IdealGasThermo(
            vib_energies,
            geometry,
            potentialenergy=self.energy,
            atoms=self.atoms,
            symmetrynumber=point_group_data.rotation_number,
            spin=spin,
            ignore_imag_modes=True,
        )
    
    def run_harmonic_thermo(self) -> HarmonicThermo:
        """
        Create a HarmonicThermo object for a molecule from a given vibrational analysis.
        This works for free gases, solids and adsorbates.

        Returns
        -------
        HarmonicThermo object
        """

        # Ensure all negative modes are made complex
        for i, f in enumerate(self.vib_freqs):
            if not isinstance(f, complex) and f < 0:
                self.vib_freqs[i] = complex(0 - f * 1j)

        # Convert vibrational frequencies to energies
        vib_energies = [f * units.invcm for f in self.vib_freqs]

        return HarmonicThermo(vib_energies=vib_energies, potentialenergy=self.energy, ignore_imag_modes=True)
