"""
Core recipes for the NewtonNet code
"""
from __future__ import annotations

from typing import Literal

import covalent as ct
import numpy as np
from ase.atoms import Atoms
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


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def relax_job(
    atoms: Atoms, newtonnet_kwargs: dict | None = None, opt_swaps: dict | None = None
) -> dict:
    """
    TODO: docstring
    """
    newtonnet_kwargs = newtonnet_kwargs or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": "sella",
        "optimizer_kwargs": {"order": 0},
    }
    opt_flags = opt_defaults | opt_swaps

    # Define calculator
    mlcalculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        config_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator

    # Run optimization
    dyn = run_ase_opt(atoms, **opt_flags)
    return summarize_opt_run(dyn, additional_fields={"name": "NewtonNet Relax"})


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def ts_job(
    atoms: Atoms,
    use_custom_hessian: bool = False,
    temperature: float = 298.15,
    pressure: float = 1.0,
    newtonnet_kwargs: dict | None = None,
    opt_swaps: dict | None = None,
) -> dict:
    """
    # TODO: docstring
    """
    newtonnet_kwargs = newtonnet_kwargs or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": "sella",
        "optimizer_kwargs": {"diag_every_n": 0} if use_custom_hessian else {},
    }
    opt_flags = opt_defaults | opt_swaps

    # Define calculator
    mlcalculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        config_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator

    if use_custom_hessian:
        if opt_flags["optimizer"].lower() != "sella":
            raise ValueError("Custom hessian can only be used with Sella.")

        atoms.calc.calculate()
        hessian = atoms.calc.results["hessian"].reshape((-1, 3 * len(atoms)))
        opt_defaults["optimizer_kwargs"]["hessian_function"] = hessian
        # TODO: I think you may need to re-initialize the calculator
        # object after this so that it's "blank" when you do
        # run_ase_opt. Please check.

    # Run the TS optimization
    dyn = run_ase_opt(atoms, **opt_flags)
    ts_summary = summarize_opt_run(dyn, additional_fields={"name": "NewtonNet TS"})

    # Run a frequency calculation
    thermo_summary = freq_job(
        ts_summary["atoms"],
        temperature=temperature,
        pressure=pressure,
        newtonnet_kwargs=newtonnet_kwargs,
    )

    return {"ts": ts_summary, "thermo": thermo_summary}


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def irc_job(
    atoms: Atoms,
    direction: Literal["forward", "reverse"] = "forward",
    temperature: float = 298.15,
    pressure: float = 1.0,
    newtonnet_kwargs: dict | None = None,
    opt_swaps: dict | None = None,
) -> dict:
    """
    TODO: docstrings
    """

    newtonnet_kwargs = newtonnet_kwargs or {}
    opt_swaps = opt_swaps or {}

    opt_defaults = {
        "fmax": 0.01,
        "max_steps": 1000,
        "optimizer": "sella_irc",
        "run_kwargs": {"direction": direction.lower()},
    }
    opt_flags = opt_defaults | opt_swaps

    # Define calculator
    mlcalculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        config_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator

    # Run IRC
    dyn = run_ase_opt(atoms, **opt_flags)
    summary_irc = summarize_opt_run(dyn, additional_fields={"name": "NewtonNet IRC"})

    # Run frequency job
    thermo_summary = freq_job(
        summary_irc["atoms"],
        temperature=temperature,
        pressure=pressure,
        newtonnet_kwargs=newtonnet_kwargs,
    )
    return {"irc": summary_irc, "thermo": thermo_summary}


@ct.electron
@requires(NewtonNet, "NewtonNet must be installed. Try pip install quacc[newtonnet]")
def quasi_irc_job(
    atoms: Atoms,
    direction: Literal["forward", "reverse"] = "forward",
    temperature: float = 298.15,
    pressure: float = 1.0,
    newtonnet_kwargs: dict | None = None,
    irc_swaps: dict | None = None,
    opt_swaps: dict | None = None,
) -> dict:
    """
    TODO: docstrings
    """
    newtonnet_kwargs = newtonnet_kwargs or {}
    irc_swaps = irc_swaps or {}
    opt_swaps = opt_swaps or {}

    irc_defaults = {
        "fmax": 0.01,
        "max_steps": 5,
        "optimizer": "sella_irc",
        "run_kwargs": {"direction": direction.lower()},
    }
    irc_flags = irc_defaults | irc_swaps

    opt_defaults = {"fmax": 0.01, "max_steps": 1000, "optimizer": "sella"}
    opt_flags = opt_defaults | opt_swaps

    # Define calculator
    mlcalculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        config_path=SETTINGS.NEWTONNET_CONFIG_PATH,
        **newtonnet_kwargs,
    )
    atoms.calc = mlcalculator

    # Run IRC
    irc_summary = irc_job(atoms, newtonnet_kwargs=newtonnet_kwargs, **irc_flags)

    # Run opt
    opt_summary = relax_job(
        irc_summary["atoms"], newtonnet_kwargs=newtonnet_kwargs, **opt_flags
    )

    # Run frequency
    thermo_summary = freq_job(
        opt_summary["atoms"],
        temperature=temperature,
        pressure=pressure,
        newtonnet_kwargs=newtonnet_kwargs,
    )

    return {"irc": irc_summary, "opt": opt_summary, "thermo": thermo_summary}


# TODO: I think it is possible to get rid of all the unit conversion functions
# and instead use the `VibrationsData` class in `ase.vibrations.data.py`
# See details here: https://gitlab.com/ase/ase/-/blob/master/ase/vibrations/data.py
# Once that class is contructed, you can do `VibrationsData.get_frequencies()`
# to return the frequencies in cm^-1, as desired. You can then get rid of
# the `_get_freq_in_cm_inv` and `_mass_weighted_hessian` functions.
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

    # TODO: If you are successful in using the `VibrationsData` class, you
    # can then return the following instead for a much richer and consistent output:
    # return {
    #     "vib": summarize_vib_run(
    #         vibrations_data, additional_fields={"name": "NewtonNet Vibrations"}
    #     ),
    #     "thermo": summarize_thermo_run(
    #         igt,
    #         temperature=temperature,
    #         pressure=pressure,
    #         additional_fields={"name": "NewtonNet Thermo"},
    #     ),
    # }

    return summarize_thermo_run(
        igt,
        temperature=temperature,
        pressure=pressure,
        additional_fields={"name": "NewtonNet Thermo"},
    )


# TODO: Can potentially replace with `VibrationsData` (see above)
def _get_freq_in_cm_inv(
    masses: np.array, reshaped_hessian: np.ndarray
) -> list[np.array, np.array]:
    # Calculate mass-weighted Hessian
    mass_weighted_hessian = _mass_weighted_hessian(masses, reshaped_hessian)

    # Calculate eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eig(mass_weighted_hessian)
    eigvals = np.sort(np.real(eigvals))

    # Calculate frequencies in cm^-1
    freqs = np.emath.sqrt(eigvals) * fs * (10**15) / (_c * 100 * 2 * np.pi)
    return freqs, eigvecs


# TODO: Can potentially replace with `VibrationsData` (see above)
def _mass_weighted_hessian(masses: np.array, hessian: np.ndarray) -> np.ndarray:
    """
    Calculates the mass-weighted Hessian matrix.

    Parameters:
        masses (array): An array of atom masses.
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

    if len(masses) != hessian.shape[0] // 3 or hessian.shape[0] != hessian.shape[1]:
        raise ValueError("Incompatible dimensions of masses and hessian.")

    sqrt_masses = np.sqrt(np.outer(masses, masses))
    return hessian / np.tile(sqrt_masses, (3, 3))
