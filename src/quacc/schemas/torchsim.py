from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

from torch_sim.runners import (
    generate_energy_convergence_fn,
    generate_force_convergence_fn,
)

if TYPE_CHECKING:
    from collections.abc import Callable


class TSModelType(StrEnum):
    """Enum for model types."""

    FAIRCHEMV1 = "FairChemV1Model"
    FAIRCHEM = "FairChemModel"
    GRAPHPESWRAPPER = "GraphPESWrapper"
    MACE = "MaceModel"
    MATTERSIM = "MatterSimModel"
    METATOMIC = "MetatomicModel"
    NEQUIPFRAMEWORK = "NequIPFrameworkModel"
    ORB = "OrbModel"
    SEVENNET = "SevenNetModel"
    LENNARD_JONES = "LennardJonesModel"


class ConvergenceFn(StrEnum):
    """Enum for convergence function types."""

    ENERGY = "energy"
    FORCE = "force"


CONVERGENCE_FN_REGISTRY: dict[str, Callable] = {
    "energy": generate_energy_convergence_fn,
    "force": generate_force_convergence_fn,
}
