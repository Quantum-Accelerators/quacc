"""Schemas for storing thermo data."""

from __future__ import annotations

from logging import getLogger
from typing import TYPE_CHECKING

from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.units import invcm
from emmet.core.symmetry import PointGroupData
from pymatgen.io.ase import AseAtomsAdaptor

from quacc import __version__, get_settings
from quacc.atoms.core import get_spin_multiplicity_attribute
from quacc.schemas.atoms import atoms_to_metadata
from quacc.utils.dicts import clean_dict
from quacc.utils.files import get_uri

if TYPE_CHECKING:
    from ase.atoms import Atoms

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import ThermoSchema

LOGGER = getLogger(__name__)


class ThermoSummarize:
    """
    Summarize an ASE thermochemistry object for storage in a database.
    """

    def __init__(
        self,
        atoms: Atoms,
        vib_freqs: list[float | complex],
        energy: float = 0.0,
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
        additional_fields
            Additional fields to store in the document.

        Returns
        -------
        None
        """
        self.atoms = atoms
        # Make sure vibrational freqs are imaginary, not negative
        vib_freqs_ = vib_freqs.copy()
        for i, f in enumerate(vib_freqs_):
            if not isinstance(f, complex) and f < 0:
                vib_freqs_[i] = complex(0 - f * 1j)
        self.vib_freqs = vib_freqs_
        self.vib_energies = [f * invcm for f in self.vib_freqs]
        self.energy = energy
        self.additional_fields = additional_fields or {}
        self._settings = get_settings()

    def ideal_gas(
        self,
        temperature: float = 298.15,
        pressure: float = 1.0,
        spin_multiplicity: int | None = None,
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
        spin_multiplicity
            Spin multiplicity of the system. If not provided, will attempt to detect
            from the Atoms object.

        Returns
        -------
        ThermoSchema
            Dictionary representation of the task document
        """

        # Get the spin multiplicity
        if spin_multiplicity is None:
            spin_multiplicity = get_spin_multiplicity_attribute(self.atoms)
            LOGGER.warning(
                f"Using a spin multiplicity of {spin_multiplicity} for IdealGasThermo."
            )

        # Generate the ASE IdealGasThermo object
        igt = self._make_ideal_gas(spin_multiplicity=spin_multiplicity)

        # Tabulate inputs and outputs
        inputs = {
            "parameters_thermo": {
                "temperature": temperature,
                "pressure": pressure,
                "sigma": igt.sigma,
                "spin_multiplicity": spin_multiplicity,
                "vib_freqs": [e / invcm for e in igt.vib_energies],
                "vib_energies": igt.vib_energies.tolist(),
                "n_imag": igt.n_imag,
                "method": "ideal_gas",
            },
            "quacc_version": __version__,
        }

        results = {
            "results": {
                "energy": igt.potentialenergy,
                "enthalpy": igt.get_enthalpy(temperature),
                "entropy": igt.get_entropy(temperature, pressure * 10**5),
                "gibbs_energy": igt.get_gibbs_energy(temperature, pressure * 10**5),
                "zpe": igt.get_ZPE_correction(),
            }
        }

        unsorted_task_doc = (
            atoms_to_metadata(igt.atoms) | inputs | results | self.additional_fields
        )
        return clean_dict(unsorted_task_doc)

    def harmonic(
        self, temperature: float = 298.15, pressure: float = 1.0
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

        Returns
        -------
        ThermoSchema
            Dictionary representation of the task document
        """
        # Generate the ASE HarmonicThermo object
        harmonic_thermo = self._make_harmonic_thermo()

        # Tabulate inputs and outputs
        inputs = {
            "parameters_thermo": {
                "temperature": temperature,
                "pressure": pressure,
                "vib_freqs": [e / invcm for e in harmonic_thermo.vib_energies],
                "vib_energies": harmonic_thermo.vib_energies.tolist(),
                "n_imag": harmonic_thermo.n_imag,
                "method": "harmonic",
            },
            "nid": get_uri(self.directory).split(":")[0],
            "dir_name": self.directory,
            "quacc_version": __version__,
        }

        results = {
            "results": {
                "energy": harmonic_thermo.potentialenergy,
                "helmholtz_energy": harmonic_thermo.get_helmholtz_energy(temperature),
                "internal_energy": harmonic_thermo.get_internal_energy(temperature),
                "entropy": harmonic_thermo.get_entropy(temperature),
                "zpe": harmonic_thermo.get_ZPE_correction(),
            }
        }

        unsorted_task_doc = (
            atoms_to_metadata(self.atoms) | inputs | results | self.additional_fields
        )
        return clean_dict(unsorted_task_doc)

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
            self.vib_energies,
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
        return HarmonicThermo(
            vib_energies=self.vib_energies,
            potentialenergy=self.energy,
            ignore_imag_modes=True,
        )
