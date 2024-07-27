"""Schemas for storing thermo data."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from ase import units
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from emmet.core.symmetry import PointGroupData
from pymatgen.io.ase import AseAtomsAdaptor

from quacc import QuaccDefault, get_settings
from quacc.atoms.core import get_spin_multiplicity_attribute
from quacc.schemas.atoms import atoms_to_metadata
from quacc.utils.dicts import finalize_dict

if TYPE_CHECKING:
    from ase.atoms import Atoms

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from maggma.core import Store

    from quacc.types import DefaultSetting, ThermoSchema

LOGGER = logging.getLogger(__name__)


class ThermoSummarize:
    """
    Summarize an ASE thermochemistry object for storage in a database.
    """

    def __init__(
        self,
        atoms: Atoms,
        vib_freqs: list[float | complex],
        energy: float = 0.0,
        directory: str | Path | None = None,
        charge_and_multiplicity: tuple[int, int] | None = None,
        additional_fields: dict[str, Any] | None = None,
    ) -> None:
        """
        Initialize the ThermoSummarize object.

        Parameters
        ----------
        atoms
            ASE Atoms object.
        vib_freqs
            Vibrational frequencies.
        energy
            Potential energy used for a reference in thermochemistry calculations.
        directory
            Directory to store the output files. Defaults to the directory of the
            atoms object's calculator, if available.
        charge_and_multiplicity
            Charge and multiplicity of the atoms object.
        additional_fields
            Additional fields to store in the document.

        Returns
        -------
        None
        """
        self.atoms = atoms
        self.vib_freqs = vib_freqs
        self.energy = energy
        self.directory = Path(directory or atoms.calc.directory)
        self.charge_and_multiplicity = charge_and_multiplicity
        self.additional_fields = additional_fields or {}
        self._settings = get_settings()

    def ideal_gas(
        self,
        temperature: float = 298.15,
        pressure: float = 1.0,
        store: Store | None | DefaultSetting = QuaccDefault,
    ) -> ThermoSchema:
        """
        Get tabulated results from an ASE IdealGasThermo object and store them in a
        database-friendly format.

        Parameters
        ----------
        temperature
            Temperature in Kelvins.
        pressure
            Pressure in bar.
        store
            Whether to store the document in the database.

        Returns
        -------
        ThermoSchema
            Dictionary representation of the task document
        """
        store = self._settings.STORE if store == QuaccDefault else store
        if self.charge_and_multiplicity:
            spin_multiplicity = self.charge_and_multiplicity[1]
        else:
            spin_multiplicity = get_spin_multiplicity_attribute(self.atoms)
            LOGGER.info(
                f"No multiplicity provided. Automatically detecting a spin multiplicity of {spin_multiplicity} from the Atoms object."
            )
        igt = self._make_ideal_gas(spin_multiplicity=spin_multiplicity)

        inputs = {
            "parameters_thermo": {
                "temperature": temperature,
                "pressure": pressure,
                "sigma": igt.sigma,
                "spin_multiplicity": spin_multiplicity,
                "vib_freqs": [e / units.invcm for e in igt.vib_energies],
                "vib_energies": igt.vib_energies.tolist(),
                "n_imag": igt.n_imag,
                "method": "ideal_gas",
            }
        }

        results = {
            "results": {
                "energy": igt.potentialenergy,
                "enthalpy": igt.get_enthalpy(temperature, verbose=self._settings.DEBUG),
                "entropy": igt.get_entropy(
                    temperature, pressure * 10**5, verbose=self._settings.DEBUG
                ),
                "gibbs_energy": igt.get_gibbs_energy(
                    temperature, pressure * 10**5, verbose=self._settings.DEBUG
                ),
                "zpe": igt.get_ZPE_correction(),
            }
        }

        atoms_metadata = atoms_to_metadata(
            igt.atoms, charge_and_multiplicity=self.charge_and_multiplicity
        )

        unsorted_task_doc = atoms_metadata | inputs | results | self.additional_fields
        return finalize_dict(
            unsorted_task_doc,
            directory=self.directory,
            gzip_file=self._settings.GZIP_FILES,
            store=store,
        )

    def harmonic(
        self,
        temperature: float = 298.15,
        pressure: float = 1.0,
        store: Store | None | DefaultSetting = QuaccDefault,
    ) -> ThermoSchema:
        """
        Get tabulated results from an ASE HarmonicThermo object and store them in a
        database-friendly format.

        Parameters
        ----------
        temperature
            Temperature in Kelvins.
        pressure
            Pressure in bar.
        store
            Whether to store the document in the database.

        Returns
        -------
        ThermoSchema
            Dictionary representation of the task document
        """
        store = self._settings.STORE if store == QuaccDefault else store
        harmonic_thermo = self._make_harmonic_thermo()

        inputs = {
            "parameters_thermo": {
                "temperature": temperature,
                "pressure": pressure,
                "vib_freqs": [e / units.invcm for e in harmonic_thermo.vib_energies],
                "vib_energies": harmonic_thermo.vib_energies.tolist(),
                "n_imag": harmonic_thermo.n_imag,
                "method": "harmonic",
            }
        }

        results = {
            "results": {
                "energy": harmonic_thermo.potentialenergy,
                "helmholtz_energy": harmonic_thermo.get_helmholtz_energy(
                    temperature, verbose=self._settings.DEBUG
                ),
                "internal_energy": harmonic_thermo.get_internal_energy(
                    temperature, verbose=self._settings.DEBUG
                ),
                "entropy": harmonic_thermo.get_entropy(
                    temperature, verbose=self._settings.DEBUG
                ),
                "zpe": harmonic_thermo.get_ZPE_correction(),
            }
        }

        atoms_metadata = atoms_to_metadata(self.atoms)
        unsorted_task_doc = atoms_metadata | inputs | results | self.additional_fields
        return finalize_dict(
            unsorted_task_doc,
            directory=self.directory,
            gzip_file=self._settings.GZIP_FILES,
            store=store,
        )

    def _make_ideal_gas(self, spin_multiplicity: int | None = None) -> IdealGasThermo:
        """
        Make an ASE IdealGasThermo object.

        Parameters
        ----------
        spin_multiplicity
            Spin multiplicity of the system. A value of None will attempt to detect
            the spin multiplicity from the Atoms object.

        Returns
        -------
        IdealGasThermo
            ASE IdealGasThermo object
        """
        # Ensure all negative modes are made complex
        for i, f in enumerate(self.vib_freqs):
            if not isinstance(f, complex) and f < 0:
                self.vib_freqs[i] = complex(0 - f * 1j)

        # Convert vibrational frequencies to energies
        vib_energies = [f * units.invcm for f in self.vib_freqs]

        # Get the spin from the Atoms object.
        spin = round((spin_multiplicity - 1) / 2, 1) if spin_multiplicity else 0

        # Get symmetry for later use
        mol = AseAtomsAdaptor().get_molecule(self.atoms, charge_spin_check=False)
        point_group_data = PointGroupData().from_molecule(mol)

        # Get the geometry
        natoms = len(self.atoms)
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

    def _make_harmonic_thermo(self) -> HarmonicThermo:
        """
        Make an ASE HarmonicThermo object.

        Returns
        -------
        HarmonicThermo
            ASE HarmonicThermo object
        """
        # Ensure all negative modes are made complex
        for i, f in enumerate(self.vib_freqs):
            if not isinstance(f, complex) and f < 0:
                self.vib_freqs[i] = complex(0 - f * 1j)

        # Convert vibrational frequencies to energies
        vib_energies = [f * units.invcm for f in self.vib_freqs]

        return HarmonicThermo(
            vib_energies=vib_energies,
            potentialenergy=self.energy,
            ignore_imag_modes=True,
        )
