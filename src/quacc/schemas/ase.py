"""Schemas for storing ASE-based data."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ase import units
from ase.io import read
from ase.vibrations.data import VibrationsData

from quacc import SETTINGS, __version__
from quacc.atoms.core import get_final_atoms_from_dyn
from quacc.schemas.atoms import atoms_to_metadata
from quacc.schemas.prep import prep_next_run
from quacc.utils.dicts import clean_task_doc, recursive_dict_merge
from quacc.utils.files import get_uri
from quacc.wflow_tools.db import results_to_db

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from ase.io import Trajectory
    from ase.optimize.optimize import Optimizer
    from ase.thermochemistry import IdealGasThermo
    from ase.vibrations import Vibrations
    from maggma.core import Store

    from quacc.schemas._aliases.ase import (
        OptSchema,
        RunSchema,
        ThermoSchema,
        VibSchema,
        VibThermoSchema,
    )


def summarize_run(
    final_atoms: Atoms,
    input_atoms: Atoms,
    charge_and_multiplicity: tuple[int, int] | None = None,
    move_magmoms: bool = False,
    additional_fields: dict[str, Any] | None = None,
    store: Store | bool | None = None,
) -> RunSchema:
    """
    Get tabulated results from an Atoms object and calculator and store them in a
    database-friendly format. This is meant to be compatible with all calculator types.

    Parameters
    ----------
    final_atoms
        ASE Atoms following a calculation. A calculator must be attached.
    input_atoms
        Input ASE Atoms object to store.
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule
        metadata.
    move_magmoms
        Whether to move the final magmoms of the original Atoms object to the
        initial magmoms of the returned Atoms object.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None,
        `SETTINGS.STORE` will be used.

    Returns
    -------
    RunSchema
        Dictionary representation of the task document
    """

    additional_fields = additional_fields or {}
    store = SETTINGS.STORE if store is None else store

    if not final_atoms.calc:
        msg = "ASE Atoms object has no attached calculator."
        raise ValueError(msg)
    if not final_atoms.calc.results:
        msg = "ASE Atoms object's calculator has no results."
        raise ValueError(msg)

    directory = final_atoms.calc.directory
    uri = get_uri(directory)
    input_atoms_metadata = (
        atoms_to_metadata(
            input_atoms,
            charge_and_multiplicity=charge_and_multiplicity,
            store_pmg=False,
        )
        if input_atoms
        else None
    )
    inputs = {
        "parameters": final_atoms.calc.parameters,
        "nid": uri.split(":")[0],
        "dir_name": directory,
        "input_atoms": input_atoms_metadata,
        "quacc_version": __version__,
    }

    results = {"results": final_atoms.calc.results}

    atoms_to_store = prep_next_run(final_atoms, move_magmoms=move_magmoms)

    if final_atoms:
        final_atoms_metadata = atoms_to_metadata(
            atoms_to_store, charge_and_multiplicity=charge_and_multiplicity
        )
    else:
        final_atoms_metadata = {}

    unsorted_task_doc = recursive_dict_merge(
        final_atoms_metadata, inputs, results, additional_fields
    )
    task_doc = clean_task_doc(unsorted_task_doc)

    if store:
        results_to_db(store, task_doc)

    return task_doc


def summarize_opt_run(
    dyn: Optimizer,
    trajectory: Trajectory | list[Atoms] = None,
    check_convergence: bool | None = None,
    charge_and_multiplicity: tuple[int, int] | None = None,
    move_magmoms: bool = False,
    additional_fields: dict[str, Any] | None = None,
    store: Store | bool | None = None,
) -> OptSchema:
    """
    Get tabulated results from an ASE Atoms trajectory and store them in a database-
    friendly format. This is meant to be compatible with all calculator types.

    Parameters
    ----------
    dyn
        ASE Optimizer object.
    trajectory
        ASE Trajectory object or list[Atoms] from reading a trajectory file. If
        None, the trajectory must be found in dyn.traj_atoms.
    check_convergence
        Whether to check the convergence of the calculation. Defaults to True in
        settings.
    charge_and_multiplicity
        Charge and spin multiplicity of the Atoms object, only used for Molecule
        metadata.
    move_magmoms
        Whether to move the final magmoms of the original Atoms object to the
        initial magmoms of the returned Atoms object.
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None,
        `SETTINGS.STORE` will be used.

    Returns
    -------
    OptSchema
        Dictionary representation of the task document
    """

    check_convergence = (
        SETTINGS.CHECK_CONVERGENCE if check_convergence is None else check_convergence
    )
    additional_fields = additional_fields or {}
    store = SETTINGS.STORE if store is None else store

    # Get trajectory
    if not trajectory:
        trajectory = (
            dyn.traj_atoms
            if hasattr(dyn, "traj_atoms")
            else read(dyn.trajectory.filename, index=":")
        )

    initial_atoms = trajectory[0]
    final_atoms = get_final_atoms_from_dyn(dyn)
    directory = final_atoms.calc.directory

    # Check convergence
    is_converged = dyn.converged()
    if check_convergence and not is_converged:
        msg = f"Optimization did not converge. Refer to {directory}"
        raise RuntimeError(msg)

    # Base task doc
    base_task_doc = summarize_run(
        final_atoms,
        initial_atoms,
        charge_and_multiplicity=charge_and_multiplicity,
        move_magmoms=move_magmoms,
        store=False,
    )

    # Clean up the opt parameters
    parameters_opt = dyn.todict()
    parameters_opt.pop("logfile", None)
    parameters_opt.pop("restart", None)

    opt_fields = {
        "fmax": getattr(dyn, "fmax", None),
        "parameters_opt": parameters_opt,
        "converged": is_converged,
        "nsteps": dyn.get_number_of_steps(),
        "trajectory": trajectory,
        "trajectory_results": [atoms.calc.results for atoms in trajectory],
    }

    # Create a dictionary of the inputs/outputs
    unsorted_task_doc = recursive_dict_merge(
        base_task_doc, opt_fields, additional_fields
    )
    task_doc = clean_task_doc(unsorted_task_doc)

    if store:
        results_to_db(store, task_doc)

    return task_doc


def summarize_vib_and_thermo(
    vib: Vibrations,
    igt: IdealGasThermo,
    temperature: float = 298.15,
    pressure: float = 1.0,
    charge_and_multiplicity: tuple[int, int] | None = None,
    additional_fields: dict[str, Any] | None = None,
    store: Store | bool | None = None,
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
        Maggma Store object to store the results in. If None,
        `SETTINGS.STORE` will be used.

    Returns
    -------
    VibThermoSchema
        A dictionary that merges the `VibSchema` and `ThermoSchema`.
    """
    additional_fields = additional_fields or {}
    store = SETTINGS.STORE if store is None else store

    vib_task_doc = _summarize_vib_run(
        vib, charge_and_multiplicity=charge_and_multiplicity, store=False
    )
    thermo_task_doc = _summarize_ideal_gas_thermo(
        igt,
        temperature=temperature,
        pressure=pressure,
        charge_and_multiplicity=charge_and_multiplicity,
        store=False,
    )

    unsorted_task_doc = recursive_dict_merge(
        vib_task_doc, thermo_task_doc, additional_fields
    )
    task_doc = clean_task_doc(unsorted_task_doc)

    if store:
        results_to_db(store, task_doc)

    return task_doc


def _summarize_vib_run(
    vib: Vibrations | VibrationsData,
    charge_and_multiplicity: tuple[int, int] | None = None,
    additional_fields: dict[str, Any] | None = None,
    store: Store | bool | None = None,
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
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None,
        `SETTINGS.STORE` will be used.

    Returns
    -------
    VibSchema
        Dictionary representation of the task document
    """
    additional_fields = additional_fields or {}
    store = SETTINGS.STORE if store is None else store

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

    unsorted_task_doc = recursive_dict_merge(
        atoms_metadata, inputs, results, additional_fields
    )
    task_doc = clean_task_doc(unsorted_task_doc)

    if store:
        results_to_db(store, task_doc)

    return task_doc


def _summarize_ideal_gas_thermo(
    igt: IdealGasThermo,
    temperature: float = 298.15,
    pressure: float = 1.0,
    charge_and_multiplicity: tuple[int, int] | None = None,
    additional_fields: dict[str, Any] | None = None,
    store: Store | bool | None = None,
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
    additional_fields
        Additional fields to add to the task document.
    store
        Maggma Store object to store the results in. If None,
        `SETTINGS.STORE` will be used.

    Returns
    -------
    ThermoSchema
        Dictionary representation of the task document
    """

    additional_fields = additional_fields or {}
    store = SETTINGS.STORE if store is None else store

    spin_multiplicity = round(2 * igt.spin + 1)

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
            "enthalpy": igt.get_enthalpy(temperature, verbose=SETTINGS.DEBUG),
            "entropy": igt.get_entropy(
                temperature, pressure * 10**5, verbose=SETTINGS.DEBUG
            ),
            "gibbs_energy": igt.get_gibbs_energy(
                temperature, pressure * 10**5, verbose=SETTINGS.DEBUG
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

    unsorted_task_doc = recursive_dict_merge(
        atoms_metadata, inputs, results, additional_fields
    )
    task_doc = clean_task_doc(unsorted_task_doc)

    if store:
        results_to_db(store, task_doc)

    return task_doc
