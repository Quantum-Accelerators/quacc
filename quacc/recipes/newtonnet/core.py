"""
Core recipes for the NewtonNet code
"""
from __future__ import annotations

import os

import covalent as ct
import numpy as np
from ase.atoms import Atoms
from ase.data import atomic_masses
from ase.units import _c, fs
from monty.dev import requires

from quacc import SETTINGS
from quacc.schemas.ase import summarize_opt_run, summarize_thermo_run
from quacc.util.calc import run_ase_opt
from quacc.util.thermo import ideal_gas

try:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
except ImportError:
    NewtonNet = None

# TODO: Refactor so there's less copy/paste
# TODO: Do we need the model_path and config_path as kwargs or are global settings okay?
# TODO: Are there "standard" models and configs provided with NewtonNet? If so,
#   let's add them to quacc.presets so the user doesn't need to download them.
# TODO: Might be better to add a convenience function to NewtonNet that returns frequencies
# in the desired cm^-1 units so we don't need to do that in Quacc
# TODO: Not sure we want all these to be individual Slurm jobs since they're fast? Might
#      be better to make them regular functions and have a single Slurm job combining them.
# TODO: Add docstrings and typehints for all functions/classes.


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def ts_job(
    atoms: Atoms,
    model_path: str = SETTINGS.NEWTONNET_MODEL_PATH,
    config_path: str = SETTINGS.NEWTONNET_CONFIG_PATH,
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
    if not model_path or not config_path:
        raise ValueError(
            "model_path and config_path must be specified in either the global Quacc settings or as kwargs."
        )

    for f in {model_path, config_path}:
        if not os.path.exists(f):
            raise ValueError(f"{f} does not exist.")

    newtonnet_kwargs = newtonnet_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    mlcalculator = NewtonNet(
        model_path=model_path,
        config_path=config_path,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator

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
        mol = summary["atoms"]
        mol.calc = mlcalculator
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

        igt = ideal_gas(
            mol,
            freqs_eV,
            energy=mlcalculator.results["energy"],
        )
        summary["thermo_results"] = summarize_thermo_run(
            igt,
            temperature=temperature,
            pressure=pressure,
            additional_fields={"name": "Sella Thermo"},
        )
    return summary


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def irc_job(
    atoms: Atoms,
    model_path: str = SETTINGS.NEWTONNET_MODEL_PATH,
    config_path: str = SETTINGS.NEWTONNET_CONFIG_PATH,
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
    for f in [model_path, config_path]:
        if not os.path.exists(f):
            raise ValueError(f"{f} does not exist.")

    newtonnet_kwargs = newtonnet_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    mlcalculator = NewtonNet(
        model_path=model_path,
        config_path=config_path,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator
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
        mol = summary["atoms"]
        mol.calc = mlcalculator
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

        igt = ideal_gas(mol, freqs_eV, energy=mlcalculator.results["energy"])
        summary["thermo_results"] = summarize_thermo_run(
            igt,
            temperature=temperature,
            pressure=pressure,
            additional_fields={"name": "Sella Thermo"},
        )
    return summary


# TODO: Must reduce number of kwargs
# TODO: Must reduce copy-pasting
@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def irc_job1(
    atoms: Atoms,
    model_path: str = SETTINGS.NEWTONNET_MODEL_PATH,
    config_path: str = SETTINGS.NEWTONNET_CONFIG_PATH,
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
    for f in [model_path, config_path]:
        if not os.path.exists(f):
            raise ValueError(f"{f} does not exist.")
    mlcalculator = NewtonNet(
        model_path=model_path,
        config_path=config_path,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator
    run_kwargs = {"direction": direction}
    dyn1 = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps1,
        optimizer=optimizer1,
        opt_kwargs=opt1_kwargs,
        run_kwargs=run_kwargs,
    )
    summary1 = summarize_opt_run(dyn1, additional_fields={"name": "Sella IRC 1"})

    atoms2 = summary1["atoms"]
    mlcalculator = NewtonNet(
        model_path=model_path,
        config_path=config_path,
        **newtonnet_kwargs,
    )
    atoms2.calc = mlcalculator
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
    summary2 = summarize_opt_run(dyn2, additional_fields={"name": "Sella IRC 2"})
    summary1["optimization"] = summary2
    if calc_final_freqs:
        mol = summary2["atoms"]
        mol.calc = mlcalculator
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

        igt = ideal_gas(mol, freqs_eV, energy=mlcalculator.results["energy"])
        summary1["thermo_results"] = summarize_thermo_run(
            igt,
            temperature=temperature,
            pressure=pressure,
            additional_fields={"name": "Sella Thermo"},
        )
    return summary1


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def freq_job(
    atoms: Atoms,
    temperature: float = 298.15,
    pressure: float = 1.0,
    newtonnet_kwargs: dict = None,
    model_path: str | list[str] = SETTINGS.NEWTONNET_MODEL_PATH,
    config_path: str | list[str] = SETTINGS.NEWTONNET_CONFIG_PATH,
) -> dict:
    """
    # TODO: Add docstrings
    """
    newtonnet_kwargs = newtonnet_kwargs or {}

    # Make sure user sets required paths
    if not model_path or not config_path:
        raise ValueError(
            "model_path and config_path must be specified in either the global Quacc settings or as kwargs."
        )

    # Define calculator
    mlcalculator = NewtonNet(
        model_path=model_path, config_path=config_path, **newtonnet_kwargs
    )
    atoms.calc = mlcalculator

    # Run calculator
    mlcalculator.calculate(atoms)
    hessian = mlcalculator.results["hessian"]

    # Calculate frequencies
    n_atoms = np.shape(hessian)[0]  # TODO: can't this be obtained from len(atoms)?
    hessian_reshaped = np.reshape(hessian, (n_atoms * 3, n_atoms * 3))
    eigvals, _ = np.linalg.eig(hessian_reshaped)
    vib_freqs = np.sqrt(eigvals)

    # TODO: Are the units right on the frequencies? cm^-1 but I don't see the conversion

    # Sort the frequencies (can remove after ASE MR 2906 is merged)
    vib_freqs = list(vib_freqs)
    vib_freqs.sort(key=np.imag, reverse=True)
    vib_freqs.sort(key=np.real)

    # Make IdealGasThermo object
    igt = ideal_gas(atoms, vib_freqs, energy=mlcalculator.results["energy"])

    return summarize_thermo_run(
        igt,
        temperature=temperature,
        pressure=pressure,
        additional_fields={"name": "Sella Thermo"},
    )


# TODO: Make docstring compatible with that used in Quacc
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
    return hessian / np.tile(sqrt_masses, (3, 3))


# TODO: type hints and docstrings
def _get_freq_in_cm_inv(masses, reshaped_hessian):
    # Calculate mass-weighted Hessian
    mass_weighted_hessian = _mass_weighted_hessian(masses, reshaped_hessian)

    # Calculate eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eig(mass_weighted_hessian)
    eigvals = np.sort(np.real(eigvals))

    # Calculate frequencies in cm^-1
    freqs = np.emath.sqrt(eigvals) * fs * (10**15) / (_c * 100 * 2 * np.pi)
    return [freqs, eigvecs]


# TODO: type hints and docstrings
def _get_hessian(atoms):
    atoms.calc.calculate(atoms)
    return atoms.calc.results["hessian"].reshape((-1, 3 * len(atoms)))
