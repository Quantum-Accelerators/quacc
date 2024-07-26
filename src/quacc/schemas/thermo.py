"""Schemas for storing thermo data."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ase import units
from ase.vibrations import Vibrations
from ase.vibrations.data import VibrationsData

from quacc import QuaccDefault, get_settings
from quacc.schemas.atoms import atoms_to_metadata
from quacc.utils.dicts import finalize_dict, recursive_dict_merge
from quacc.utils.files import get_uri

if TYPE_CHECKING:
    from typing import Any

    from ase.thermochemistry import IdealGasThermo
    from maggma.core import Store

    from quacc.types import ThermoSchema, VibSchema, VibThermoSchema


def summarize_vib_and_thermo(
    vib: Vibrations,
    igt: IdealGasThermo,
    temperature: float = 298.15,
    pressure: float = 1.0,
    charge_and_multiplicity: tuple[int, int] | None = None,
    additional_fields: dict[str, Any] | None = None,
    store: Store | None = QuaccDefault,
) -> VibThermoSchema:
    """
    Get tabulated results from an ASE Vibrations run and ASE IdealGasThermo object and
    store them in a database-friendly format.

    Parameters
    ----------
    vib
        ASE Vibrations object.
    igt
        ASE IdealGasThermo object.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule
        metadata.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. Defaults to  `QuaccSettings.STORE`.

    Returns
    -------
    VibThermoSchema
        A dictionary that merges the `VibSchema` and `ThermoSchema`.
    """
    settings = get_settings()
    store = settings.STORE if store == QuaccDefault else store

    vib_task_doc = _summarize_vib_run(
        vib, charge_and_multiplicity=charge_and_multiplicity
    )
    thermo_task_doc = _summarize_ideal_gas_thermo(
        igt,
        temperature=temperature,
        pressure=pressure,
        charge_and_multiplicity=charge_and_multiplicity,
    )

    unsorted_task_doc = recursive_dict_merge(
        vib_task_doc, thermo_task_doc, additional_fields
    )

    return finalize_dict(
        unsorted_task_doc,
        vib.atoms.calc.directory if isinstance(vib, Vibrations) else None,
        gzip_file=settings.GZIP_FILES,
        store=store,
    )


def _summarize_vib_run(
    vib: Vibrations | VibrationsData,
    charge_and_multiplicity: tuple[int, int] | None = None,
) -> VibSchema:
    """
    Get tabulated results from an ASE Vibrations object and store them in a database-
    friendly format.

    Parameters
    ----------
    vib
        ASE Vibrations object.
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule
        metadata.

    Returns
    -------
    VibSchema
        Dictionary representation of the task document
    """
    vib_freqs_raw = vib.get_frequencies().tolist()
    vib_energies_raw = vib.get_energies().tolist()

    # Convert imaginary modes to negative values for DB storage
    for i, f in enumerate(vib_freqs_raw):
        if np.imag(f) > 0:
            vib_freqs_raw[i] = -np.abs(f)
            vib_energies_raw[i] = -np.abs(vib_energies_raw[i])
        else:
            vib_freqs_raw[i] = np.abs(f)
            vib_energies_raw[i] = np.abs(vib_energies_raw[i])

    if isinstance(vib, VibrationsData):
        atoms = vib._atoms
        inputs = {}
    else:
        atoms = vib.atoms
        directory = atoms.calc.directory
        uri = get_uri(directory)

        inputs = {
            "parameters": atoms.calc.parameters,
            "parameters_vib": {
                "delta": vib.delta,
                "direction": vib.direction,
                "method": vib.method,
                "ndof": vib.ndof,
                "nfree": vib.nfree,
            },
            "nid": uri.split(":")[0],
            "dir_name": directory,
        }

    atoms_metadata = atoms_to_metadata(
        atoms, charge_and_multiplicity=charge_and_multiplicity
    )

    # Get the true vibrational modes
    natoms = len(atoms)
    if natoms == 1:
        vib_freqs = []
        vib_energies = []
    elif atoms.pbc.any():
        vib_freqs = vib_freqs_raw
        vib_energies = vib_energies_raw
    else:
        # Sort by absolute value
        vib_freqs_raw_sorted = vib_freqs_raw.copy()
        vib_energies_raw_sorted = vib_energies_raw.copy()
        vib_freqs_raw_sorted.sort(key=np.abs)
        vib_energies_raw_sorted.sort(key=np.abs)

        # Cut the 3N-5 or 3N-6 modes based on their absolute value
        n_modes = (
            3 * natoms - 5 if atoms_metadata["symmetry"]["linear"] else 3 * natoms - 6
        )
        vib_freqs = vib_freqs_raw_sorted[-n_modes:]
        vib_energies = vib_energies_raw_sorted[-n_modes:]

    imag_vib_freqs = [f for f in vib_freqs if f < 0]

    results = {
        "results": {
            "imag_vib_freqs": imag_vib_freqs,
            "n_imag": len(imag_vib_freqs),
            "vib_energies": vib_energies,
            "vib_freqs": vib_freqs,
            "vib_energies_raw": vib_energies_raw,
            "vib_freqs_raw": vib_freqs_raw,
        }
    }

    return atoms_metadata | inputs | results


def _summarize_ideal_gas_thermo(
    igt: IdealGasThermo,
    temperature: float = 298.15,
    pressure: float = 1.0,
    charge_and_multiplicity: tuple[int, int] | None = None,
) -> ThermoSchema:
    """
    Get tabulated results from an ASE IdealGasThermo object and store them in a
    database-friendly format.

    Parameters
    ----------
    igt
        ASE IdealGasThermo object.
    temperature
        Temperature in Kelvins.
    pressure
        Pressure in bar.
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule
        metadata.

    Returns
    -------
    ThermoSchema
        Dictionary representation of the task document
    """
    spin_multiplicity = round(2 * igt.spin + 1)
    settings = get_settings()

    inputs = {
        "parameters_thermo": {
            "temperature": temperature,
            "pressure": pressure,
            "sigma": igt.sigma,
            "spin_multiplicity": spin_multiplicity,
            "vib_freqs": [e / units.invcm for e in igt.vib_energies],
            "vib_energies": igt.vib_energies.tolist(),
            "n_imag": igt.n_imag,
        }
    }

    results = {
        "results": {
            "energy": igt.potentialenergy,
            "enthalpy": igt.get_enthalpy(temperature, verbose=settings.DEBUG),
            "entropy": igt.get_entropy(
                temperature, pressure * 10**5, verbose=settings.DEBUG
            ),
            "gibbs_energy": igt.get_gibbs_energy(
                temperature, pressure * 10**5, verbose=settings.DEBUG
            ),
            "zpe": igt.get_ZPE_correction(),
        }
    }

    if charge_and_multiplicity and spin_multiplicity != charge_and_multiplicity[1]:
        msg = (
            "The IdealGasThermo spin multiplicity does not match the user-specified multiplicity.",
        )
        raise ValueError(msg)

    atoms_metadata = atoms_to_metadata(
        igt.atoms, charge_and_multiplicity=charge_and_multiplicity
    )

    return atoms_metadata | inputs | results
