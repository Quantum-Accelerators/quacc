"""Units conversion functions for Molecular Dynamics"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
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

    if "friction" in converted_kwargs:
        converted_kwargs["friction"] = converted_kwargs.pop("friction") / fs

    if "externalstress" in converted_kwargs:
        converted_kwargs["externalstress"] = (
            np.array(converted_kwargs.pop("externalstress")) * GPa
        )

    if "timestep" in converted_kwargs:
        converted_kwargs["timestep"] = converted_kwargs.pop("timestep") * fs

    if "taut" in converted_kwargs:
        converted_kwargs["taut"] = converted_kwargs.pop("taut") * fs

    if "taup" in converted_kwargs:
        converted_kwargs["taup"] = converted_kwargs.pop("taup") * fs

    # Pressure units
    if "pfactor" in converted_kwargs:
        converted_kwargs["pfactor"] = converted_kwargs.pop("pfactor") * GPa

    # Compressibility units
    if "compressibility_au" in converted_kwargs:
        converted_kwargs["compressibility"] = (
            converted_kwargs.pop("compressibility_au") / GPa
        )

    return converted_kwargs
