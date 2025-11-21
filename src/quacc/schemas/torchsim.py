from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

import torch_sim as ts

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
    "energy": ts.generate_energy_convergence_fn,
    "force": ts.generate_force_convergence_fn,
}


class PropertyFn(StrEnum):
    """
    Because we are not able to pass live python functions through workflow serialization,
    it is necessary to have an alternative mechanism. While the functions included here are
    quite basic, this gives users a place to patch in their own functions while maintaining
    compatibility.
    """

    ENERGY = "energy"
    FORCES = "forces"
    STRESS = "stress"
    KINETIC_ENERGY = "kinetic_energy"
    TEMPERATURE = "temperature"


PROPERTY_FN_REGISTRY: dict[str, Callable] = {
    "potential_energy": lambda state: state.energy,
    "forces": lambda state: state.forces,
    "stress": lambda state: state.stress,
    "kinetic_energy": lambda state: ts.calc_kinetic_energy(
        velocities=state.velocities, masses=state.masses
    ),
    "temperature": lambda state: state.calc_temperature(),
}
