"""Units conversion functions for Molecular Dynamics"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from ase.units import GPa, fs

if TYPE_CHECKING:
    from typing import Any

LOGGER = logging.getLogger(__name__)


def md_units(dynamics_kwargs: dict[str, Any]) -> dict[str, Any]:
    """
    Convert units for molecular dynamics parameters to ensure consistency.
    Quacc ALWAYS uses the following units:
    - Time: femtoseconds (fs)
    - Pressure: GPa
    - Temperature: Kelvin (K)
    - Compressibility: 1/GPa

    Parameters
    ----------
    dynamics_kwargs: dict[str, Any]
        Dictionary of keyword arguments for the molecular dynamics calculation.

    Returns
    -------
    dict[str, Any]
        Dictionary with unit-converted keyword arguments.
    """
    converted_kwargs = dynamics_kwargs.copy()

    # Time units
    if "ttime" in converted_kwargs:
        converted_kwargs["ttime"] = converted_kwargs.pop("ttime") * fs
        LOGGER.info(f"Converted 'ttime' to {converted_kwargs['ttime']} fs")

    if "friction" in converted_kwargs:
        converted_kwargs["friction"] = converted_kwargs.pop("friction") / fs
        LOGGER.info(f"Converted 'friction' to {converted_kwargs['friction']} fs^-1")

    if "taut" in converted_kwargs:
        converted_kwargs["taut"] = converted_kwargs.pop("taut") * fs
        LOGGER.info(f"Converted 'taut' to {converted_kwargs['taut']} fs")

    if "taup" in converted_kwargs:
        converted_kwargs["taup"] = converted_kwargs.pop("taup") * fs
        LOGGER.info(f"Converted 'taup' to {converted_kwargs['taup']} fs")

    # Pressure units
    if "pfactor" in converted_kwargs:
        converted_kwargs["pfactor"] = converted_kwargs.pop("pfactor") * GPa
        LOGGER.info(f"Converted 'pfactor' to {converted_kwargs['pfactor']} GPa")

    # Compressibility units
    if "compressibility_au" in converted_kwargs:
        converted_kwargs["compressibility"] = (
            converted_kwargs.pop("compressibility_au") / GPa
        )
        LOGGER.info(
            f"Converted 'compressibility_au' to {converted_kwargs['compressibility']} GPa^-1"
        )

    return converted_kwargs
