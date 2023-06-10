"""
Core recipes for the NewtonNet code
"""
from __future__ import annotations

import os

import covalent as ct
import numpy as np
from ase.atoms import Atoms
from ase.data import atomic_masses
from ase.units import fs
from monty.dev import requires
from scipy.constants import c

from quacc.schemas.ase import summarize_opt_run
from quacc.util.calc import run_ase_opt
from quacc.util.thermo import ideal_gas

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None

# TODO: Please all units based on ASE for internal consistency
# TODO: Add `model_path` and `settings_path` to the global Quacc settings. Then remove them from the kwargs.
# TODO: Add docstrings and typehints for all functions/classes.


@ct.electron
@requires(
    NewtonNet,
    "newtonnet-python must be installed. "
    "Check out https://github.com/ericyuan00000/NewtonNet/tree/main",
)
def ts_job(
    atoms: Atoms,
    model_path: str = "training_18/models/best_model_state.tar",
    settings_path: str = "training_18/run_scripts/config0.yml",
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "sella",
    temperature: float = 298.15,
    pressure: float = 1.0,
    ts_type: int = 1,
    calc_final_freqs: bool = True,
    newtonnet_kwargs: dict | None = None,
    opt_kwargs: dict | None = None,
) -> dict:
    for f in [model_path, settings_path]:
        if not os.path.exists(f):
            raise ValueError(f"{f} does not exist.")

    newtonnet_kwargs = newtonnet_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    atoms.calc = NewtonNet(
        model_path=model_path,
        settings_path=settings_path,
        **newtonnet_kwargs,
    )

    opt_kwargs["internal"] = True
    if ts_type == 0:
        opt_kwargs["hessian_function"] = _get_hessian(atoms)
        opt_kwargs["diag_every_n"] = 0

    dyn = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
    )
    summary = summarize_opt_run(dyn, additional_fields={"name": "Sella TS"})
    if calc_final_freqs:
        mol = traj[-1]
        mlcalculator.calculate(mol)
        hessian = mlcalculator.results["hessian"]
        n_atoms = np.shape(hessian)[0]
        hessian_reshaped = np.reshape(hessian, (n_atoms * 3, n_atoms * 3))

        atomic_numbers = mol.get_atomic_numbers()
        masses = [atomic_masses[number] for number in atomic_numbers]

        freqs_cm_inv, eigvecs = _get_freq_in_cm_inv(masses, hessian_reshaped)

        summary["frequencies_real_cm_inv"] = np.real(freqs_cm_inv)
        summary["frequencies_imag_cm_inv"] = np.imag(freqs_cm_inv)
        summary["normal_modes_real"] = np.real(eigvecs)
        summary["normal_modes_imag"] = np.imag(eigvecs)

        freqs_eV = freqs_cm_inv * 1.23981 * 10 ** (-4)
        summary["frequencies_real_eV"] = np.real(freqs_eV)
        summary["frequencies_imag_eV"] = np.imag(freqs_eV)

        thermo_summary = ideal_gas(
            mol,
            freqs_eV,
            temperature=temperature,
            pressure=pressure,
            energy=mlcalculator.results["energy"],
        )
        summary["thermo_results"] = thermo_summary["results"]
    return summary


@ct.electron
@requires(
    NewtonNet,
    "newtonnet-python must be installed. Check out "
    "https://github.com/ericyuan00000/NewtonNet/tree/main",
)
def irc_job(
    atoms: Atoms,
    model_path: str = "training_18/models/best_model_state.tar",
    settings_path: str = "training_18/run_scripts/config0.yml",
    direction: str = "forward",
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "sella_irc",
    temperature: float = 298.15,
    pressure: float = 1.0,
    calc_final_freqs: bool = True,
    newtonnet_kwargs: dict | None = None,
    opt_kwargs: dict | None = None,
) -> dict:
    for f in [model_path, settings_path]:
        if not os.path.exists(f):
            raise ValueError(f"{f} does not exist.")

    newtonnet_kwargs = newtonnet_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    mlcalculator = NewtonNet(
        model_path=model_path,
        settings_path=settings_path,
        **newtonnet_kwargs,
    )
    atoms.set_calculator(mlcalculator)
    run_kwargs = {"direction": direction}
    dyn = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
        run_kwargs=run_kwargs,
    )
    summary = summarize_opt_run(dyn, additional_fields={"name": "Sella IRC"})
    if calc_final_freqs:
        mol = traj[-1]
        mlcalculator.calculate(mol)
        hessian = mlcalculator.results["hessian"]
        n_atoms = np.shape(hessian)[0]
        hessian_reshaped = np.reshape(hessian, (n_atoms * 3, n_atoms * 3))

        atomic_numbers = mol.get_atomic_numbers()
        masses = [atomic_masses[number] for number in atomic_numbers]

        freqs_cm_inv, eigvecs = _get_freq_in_cm_inv(masses, hessian_reshaped)

        summary["frequencies_real_cm_inv"] = np.real(freqs_cm_inv)
        summary["frequencies_imag_cm_inv"] = np.imag(freqs_cm_inv)
        summary["normal_modes_real"] = np.real(eigvecs)
        summary["normal_modes_imag"] = np.imag(eigvecs)

        freqs_eV = freqs_cm_inv * 1.23981 * 10 ** (-4)
        summary["frequencies_real_eV"] = np.real(freqs_eV)
        summary["frequencies_imag_eV"] = np.imag(freqs_eV)

        thermo_summary = ideal_gas(
            mol,
            freqs_eV,
            temperature=temperature,
            pressure=pressure,
            energy=mlcalculator.results["energy"],
        )
        summary["thermo_results"] = thermo_summary["results"]
    return summary


@ct.electron
@requires(
    NewtonNet,
    "newtonnet-python must be installed. Check out "
    "https://github.com/ericyuan00000/NewtonNet/tree/main",
)
def irc_job1(
    atoms: Atoms,
    model_path: str = "training_18/models/best_model_state.tar",
    settings_path: str = "training_18/run_scripts/config0.yml",
    direction: str = "forward",
    fmax: float = 0.01,
    max_steps1: int = 5,
    max_steps2: int = 1000,
    optimizer1: str = "sella_irc",
    optimizer2: str = "sella",
    temperature: float = 298.15,
    pressure: float = 1.0,
    calc_final_freqs: bool = True,
    newtonnet_kwargs: dict | None = None,
    opt1_kwargs: dict | None = None,
    opt2_kwargs: dict | None = None,
) -> dict:
    for f in [model_path, settings_path]:
        if not os.path.exists(f):
            raise ValueError(f"{f} does not exist.")
    mlcalculator = NewtonNet(
        model_path=model_path,
        settings_path=settings_path,
        **newtonnet_kwargs,
    )
    atoms.set_calculator(mlcalculator)
    run_kwargs = {"direction": direction}
    dyn1 = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps1,
        optimizer=optimizer1,
        opt_kwargs=opt1_kwargs,
        run_kwargs=run_kwargs,
    )
    summary1 = summarize_opt_run(dyn1, additional_fields={"name": name})

    atoms2 = traj1[-1]
    mlcalculator = NewtonNet(
        model_path=model_path,
        settings_path=settings_path,
        **newtonnet_kwargs,
    )
    atoms2.set_calculator(mlcalculator)
    opt2_kwargs["internal"] = True
    opt2_kwargs["order"] = 0
    opt2_kwargs["trajectory"] = "opt2.traj"
    opt2_kwargs["restart"] = "opt2.pckl"

    dyn2 = run_ase_opt(
        atoms2,
        fmax=fmax,
        max_steps=max_steps2,
        optimizer=optimizer2,
        opt_kwargs=opt2_kwargs,
    )
    summary2 = summarize_opt_run(dyn2, additional_fields={"name": name})
    summary1["optimization"] = summary2
    if calc_final_freqs:
        mol = traj2[-1]
        mlcalculator.calculate(mol)
        hessian = mlcalculator.results["hessian"]
        n_atoms = np.shape(hessian)[0]
        hessian_reshaped = np.reshape(hessian, (n_atoms * 3, n_atoms * 3))

        atomic_numbers = mol.get_atomic_numbers()
        masses = [atomic_masses[number] for number in atomic_numbers]

        freqs_cm_inv, eigvecs = _get_freq_in_cm_inv(masses, hessian_reshaped)

        summary1["optimization"]["frequencies_real_cm_inv"] = np.real(freqs_cm_inv)
        summary1["optimization"]["frequencies_imag_cm_inv"] = np.imag(freqs_cm_inv)
        summary1["optimization"]["normal_modes_real"] = np.real(eigvecs)
        summary1["optimization"]["normal_modes_imag"] = np.imag(eigvecs)

        freqs_eV = freqs_cm_inv * 1.23981 * 10 ** (-4)
        summary1["optimization"]["frequencies_real_eV"] = np.real(freqs_eV)
        summary1["optimization"]["frequencies_imag_eV"] = np.imag(freqs_eV)

        thermo_summary = ideal_gas(
            mol,
            freqs_eV,
            temperature=temperature,
            pressure=pressure,
            energy=mlcalculator.results["energy"],
        )
        summary1["thermo_results"] = thermo_summary["results"]
    return summary1


@ct.electron
@requires(
    NewtonNet,
    "newtonnet-python must be installed. Check out "
    "https://github.com/ericyuan00000/NewtonNet/tree/main",
)
def freq_job(
    atoms: Atoms,
    model_path: str = "training_18/models/best_model_state.tar",
    settings_path: str = "training_18/run_scripts/config0.yml",
    temperature: float = 298.15,
    pressure: float = 1.0,
    calc_final_freqs: bool = True,
    newtonnet_kwargs: dict = None,
) -> dict:
    for f in [model_path, settings_path]:
        if not os.path.exists(f):
            raise ValueError(f"{f} does not exist.")

    mlcalculator = NewtonNet(
        model_path=model_path,
        settings_path=settings_path,
        **newtonnet_kwargs,
    )
    atoms.set_calculator(mlcalculator)
    mlcalculator.calculate(atoms)
    hessian = mlcalculator.results["hessian"]
    n_atoms = np.shape(hessian)[0]
    hessian_reshaped = np.reshape(hessian, (n_atoms * 3, n_atoms * 3))
    eigvals, eigvecs = np.linalg.eig(hessian_reshaped)
    # Note: frequencies are added below as np.sqrt(eigvals)
    thermo_summary = ideal_gas(
        atoms,
        np.sqrt(eigvals),
        temperature=temperature,
        pressure=pressure,
        energy=mlcalculator.results["energy"],
    )
    thermo_summary["name"] = "Sella Vibrations"
    return thermo_summary


def _mass_weighted_hessian(masses: list | np.array, hessian: np.ndarray) -> np.ndarray:
    """
    Calculates the mass-weighted Hessian matrix.

    Parameters:
        masses (array-like): A list or array of atom masses.
        hessian (ndarray): A 2D numpy array representing the Hessian matrix.

    Returns:
        mass_weighted_hessian (ndarray): A 2D numpy array representing the mass-weighted Hessian matrix.

    Raises:
        ValueError: If the dimensions of masses and hessian are not compatible.

    The mass-weighted Hessian matrix is calculated by dividing each element of the Hessian matrix
    by the square root of the product of the masses of the corresponding atoms. The resulting matrix
    represents the second derivatives of the potential energy surface, taking into account the masses
    of the atoms.

    The input hessian is assumed to be a square matrix of size 3N x 3N, where N is the number of atoms.
    The input masses should have N elements, where each element corresponds to the mass of an atom.

    The returned mass_weighted_hessian matrix will have the same dimensions as the input hessian.

    Example:
        masses = [12.01, 1.008, 1.008, 16.00, 1.008, 1.008, 1.008]
        hessian = np.zeros((21, 21))  # Example Hessian matrix
        mass_weighted_hessian = _mass_weighted_hessian(masses, hessian)
    """
    masses = np.asarray(masses)
    hessian = np.asarray(hessian)

    if len(masses) != hessian.shape[0] // 3 or hessian.shape[0] != hessian.shape[1]:
        raise ValueError("Incompatible dimensions of masses and hessian.")

    sqrt_masses = np.sqrt(np.outer(masses, masses))
    mass_weighted_hessian = hessian / np.tile(sqrt_masses, (3, 3))

    return mass_weighted_hessian


def _get_freq_in_cm_inv(masses, reshaped_hessian):
    # Calculate mass-weighted Hessian
    mass_weighted_hessian = _mass_weighted_hessian(masses, reshaped_hessian)

    # Calculate eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eig(mass_weighted_hessian)
    eigvals = np.sort(np.real(eigvals))

    # Calculate frequencies in cm^-1
    freqs = np.emath.sqrt(eigvals) * fs * (10**15) / (c * 100 * 2 * np.pi)
    return [freqs, eigvecs]


def _get_hessian(atoms):
    atoms.calc.calculate(atoms)
    return atoms.calc.results["hessian"].reshape((-1, 3 * len(atoms)))
