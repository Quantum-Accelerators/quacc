from __future__ import annotations

from typing import Any, TypedDict

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


class DosProperties(TypedDict):
    """Type hint associated with PhononSchema."""

    frequency_points: NDArray
    total_dos: NDArray


class PhononResults(TypedDict):
    thermal_properties: ThermalProperties
    mesh_properties: MeshProperties
    total_dos: DosProperties
    force_constants: NDArray


class PhonopyMetadata(TypedDict):
    """Type hint associated with PhononSchema."""

    version: str


class PhononSchema(AtomsSchema):
    """Type hint associated with [quacc.schemas.phonons.summarize_phonopy][]"""

    parameters: dict[str, Any] | None
    nid: str
    dir_name: str
    phonopy_metadata: PhonopyMetadata
    results: PhononResults
    quacc_version: str
