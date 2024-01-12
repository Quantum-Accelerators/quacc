from __future__ import annotations

from typing import TypedDict

from numpy.typing import NDArray

from quacc.schemas._aliases.ase import AtomsSchema


class ThermalProperties(TypedDict):
    """Type hint associated with PhononSchema."""

    temperatures: NDArray
    free_energy: NDArray
    entropy: NDArray
    heat_capacity: NDArray


class MeshProperties(TypedDict):
    """Type hint associated with PhononSchema."""

    qpoints: NDArray
    weights: NDArray
    frequencies: NDArray
    eigenvectors: NDArray
    group_velocities: NDArray


class PhononResults(TypedDict):
    thermal_properties: ThermalProperties
    mesh_properties: MeshProperties
    force_constants: NDArray


class PhononSchema(AtomsSchema):
    """Type hint associated with `quacc.schemas.phonons.summarize_phonopy`"""

    results: PhononResults
