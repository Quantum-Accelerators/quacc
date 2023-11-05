from __future__ import annotations

from typing import TypedDict

from numpy.typing import NDArray
from phonopy import Phonopy


class ThermalProperties(TypedDict):
    """Type hint associated with PhononSchema"""

    temperatures: NDArray
    free_energy: NDArray
    entropy: NDArray
    heat_capacity: NDArray


class PhononSchema(TypedDict):
    """Type hint associated with `quacc.schemas.phonopy.summarize_phonopy`"""

    phonon: Phonopy
    thermal_properties: ThermalProperties
