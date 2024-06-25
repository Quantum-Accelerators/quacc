"""Units conversion functions for Molecular Dynamics"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from ase.units import GPa, fs

if TYPE_CHECKING:
    from typing import Any

QUACC_BASE_UNITS = {
    "ttime": fs,
    "pressure_au": GPa,
    "temperature_K": 1.0,
    "compressibility_au": 1 / GPa,
    "timestep": fs,
    "taut": fs,
    "taup": fs,
    "pfactor": GPa,
    "externalstress": GPa,
    "friction": 1 / fs,
}


def convert_md_units(
    dynamics_kwargs: dict[str, Any], inverse: bool = False
) -> dict[str, Any]:
    """
    Convert units for molecular dynamics parameters to ensure consistency.

    Quacc always uses the following units:

    - Time: femtoseconds (fs)
    - Pressure: GPa
    - Temperature: Kelvin (K)
    - Compressibility: 1/GPa

    Parameters
    ----------
    dynamics_kwargs
        Dictionary of keyword arguments for the molecular dynamics calculation.
    inverse
        If True, convert from Quacc units to ASE units.

    Returns
    -------
    dict[str, Any]
        Dictionary with unit-converted keyword arguments.
    """
    converted_kwargs = dynamics_kwargs.copy()

    for key in converted_kwargs:
        if key in QUACC_BASE_UNITS:
            if inverse:
                converted_kwargs[key] = (
                    np.array(converted_kwargs[key]) / QUACC_BASE_UNITS[key]
                )

            else:
                converted_kwargs[key] = (
                    np.array(converted_kwargs[key]) * QUACC_BASE_UNITS[key]
                )
    return converted_kwargs
