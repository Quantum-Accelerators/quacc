"""
Core recipes for the NewtonNet code
"""
from __future__ import annotations

import os
from typing import Literal

import covalent as ct
import numpy as np
from ase.atoms import Atoms
from ase.data import atomic_masses
from ase.optimize.optimize import Optimizer
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


# TODO: Add reasonable defaults to settings
# TODO: Reduce number of kwargs
# TODO: Make sure the calculation fails gracefully if paths aren't set correctly.
# TODO: Not sure we want all these to be individual Slurm jobs since they're fast? Might
#      be better to make them regular functions and have a single Slurm job combining them.
# TODO: Add docstrings and typehints for all functions/classes.


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def ts_job(
    atoms: Atoms,
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = "sella",
    temperature: float = 298.15,
    pressure: float = 1.0,
    custom_hessian: bool = False,
    run_freq: bool = True,
    newtonnet_kwargs: dict | None = None,
    opt_kwargs: dict | None = None,
) -> dict:
    """
    # TODO: docstring
    """
    newtonnet_kwargs = newtonnet_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    # Define calculator
    mlcalculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        config_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator

    # Sella-related parameters
    if optimizer.lower() == "sella":
        opt_kwargs["internal"] = True
        if custom_hessian:
            opt_kwargs["diag_every_n"] = 0
            opt_kwargs["hessian_function"] = _get_hessian(atoms)
    elif custom_hessian:
        raise ValueError("Custom hessian can only be used with Sella.")

    # Run the optimization
    dyn = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
    )
    summary = summarize_opt_run(dyn, additional_fields={"name": "Sella TS"})

    # Run an optional frequency calculation
    if run_freq:
        freq_summary = freq_job(
            summary["atoms"],
            temperature=temperature,
            pressure=pressure,
            newtonnet_kwargs=newtonnet_kwargs,
        )
        summary["thermo_results"] = freq_summary
    return summary


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def irc_job(
    atoms: Atoms,
    direction: str = "forward",
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: str = Literal["sella_irc"],
    temperature: float = 298.15,
    pressure: float = 1.0,
    run_freq: bool = True,
    newtonnet_kwargs: dict | None = None,
    opt_kwargs: dict | None = None,
) -> dict:
    """
    TODO: docstrings
    """

    newtonnet_kwargs = newtonnet_kwargs or {}
    opt_kwargs = opt_kwargs or {}

    # Define calculator
    mlcalculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        config_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator

    # Run IRC
    run_kwargs = {"direction": direction.lower()}
    dyn = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps,
        optimizer=optimizer,
        opt_kwargs=opt_kwargs,
        run_kwargs=run_kwargs,
    )
    summary = summarize_opt_run(dyn, additional_fields={"name": "Sella IRC"})

    # Run optional frequency job
    if run_freq:
        freq_summary = freq_job(
            summary["atoms"],
            temperature=temperature,
            pressure=pressure,
            newtonnet_kwargs=newtonnet_kwargs,
        )
        summary["thermo_results"] = freq_summary
    return summary


# TODO: Must reduce number of kwargs
# TODO: Must reduce copy-pasting
@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def quasi_irc_job(
    atoms: Atoms,
    direction: str = "forward",
    fmax: float = 0.01,
    max_steps1: int = 5,
    max_steps2: int = 1000,
    optimizer1: str = Literal["sella_irc"],
    optimizer2: str = "sella",
    temperature: float = 298.15,
    pressure: float = 1.0,
    run_freq: bool = True,
    newtonnet_kwargs: dict | None = None,
    opt1_kwargs: dict | None = None,
    opt2_kwargs: dict | None = None,
) -> dict:
    """
    TODO: docstrings
    """
    newtonnet_kwargs = newtonnet_kwargs or {}
    opt1_kwargs = opt1_kwargs or {}
    opt2_kwargs = opt2_kwargs or {}

    # Define calculator
    mlcalculator1 = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        config_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator1

    # Run IRC
    run_kwargs = {"direction": direction.lower()}
    dyn = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps1,
        optimizer=optimizer1,
        opt_kwargs=opt1_kwargs,
        run_kwargs=run_kwargs,
    )
    irc_summary = summarize_opt_run(dyn, additional_fields={"name": "Sella IRC"})

    # Run opt
    atoms = irc_summary["atoms"]
    mlcalculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        config_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator
    opt2_kwargs["internal"] = True
    opt2_kwargs["order"] = 0

    dyn2 = run_ase_opt(
        atoms,
        fmax=fmax,
        max_steps=max_steps2,
        optimizer=optimizer2,
        opt_kwargs=opt2_kwargs,
    )
    opt_summary = summarize_opt_run(dyn2, additional_fields={"name": "Sella Opt"})

    summary = {"IRC": irc_summary, "Opt": opt_summary}

    # Run optional frequency
    if run_freq:
        freq_summary = freq_job(
            opt_summary["atoms"],
            temperature=temperature,
            pressure=pressure,
            newtonnet_kwargs=newtonnet_kwargs,
        )
        summary["freq"] = freq_summary

    return summary


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def freq_job(
    atoms: Atoms,
    temperature: float = 298.15,
    pressure: float = 1.0,
    newtonnet_kwargs: dict = None,
) -> dict:
    """
    # TODO: Add docstrings
    """
    newtonnet_kwargs = newtonnet_kwargs or {}

    # Define calculator
    mlcalculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        config_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator


    # Run calculator
    mlcalculator.calculate(atoms)
    hessian = mlcalculator.results["hessian"]

    # Calculate frequencies
    n_atoms = len(atoms)
    hessian_reshaped = np.reshape(hessian, (n_atoms * 3, n_atoms * 3))
    freqs_cm_inv, _ = _get_freq_in_cm_inv(atoms.get_masses(), hessian_reshaped)

    # Sort the frequencies (can remove after ASE MR 2906 is merged)
    freqs_cm_inv = list(freqs_cm_inv)
    freqs_cm_inv.sort(key=np.imag, reverse=True)
    freqs_cm_inv.sort(key=np.real)

    # Make IdealGasThermo object
    igt = ideal_gas(atoms, freqs_cm_inv, energy=mlcalculator.results["energy"])

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
    atoms.calc.calculate()
    return atoms.calc.results["hessian"].reshape((-1, 3 * len(atoms)))
