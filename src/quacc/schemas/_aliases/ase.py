"""Aliases for type hinting `quacc.schemas.ase`"""

from __future__ import annotations

from typing import Any, TypedDict

from ase.atoms import Atoms
from numpy.typing import NDArray

from quacc.schemas._aliases.atoms import AtomsSchema

Results = dict[str, Any]  # from atoms.calc.results
Parameters = dict[str, Any]  # from atoms.calc.parameters


class RunSchema(AtomsSchema):
    """Schema for [quacc.schemas.ase.summarize_run][]"""

    input_atoms: AtomsSchema | None
    nid: str
    dir_name: str
    parameters: Parameters
    results: Results
    quacc_version: str


class OptSchema(RunSchema):
    """Schema for [quacc.schemas.ase.summarize_opt_run][]"""

    fmax: float | None
    parameters_opt: dict[str, Any]  # from Optimizer.todict()
    converged: bool
    nsteps: int
    trajectory: list[Atoms]
    trajectory_results: list[Results]


class ParametersVib(TypedDict):
    delta: float
    direction: str
    method: str
    ndof: int
    nfree: int


class VibResults(TypedDict):
    imag_vib_freqs: int
    n_imag: int
    vib_energies: list[float]
    vib_freqs: list[float]
    vib_energies_raw: list[float]
    vib_freqs_raw: list[float]


class VibSchema(AtomsSchema):
    parameters: Parameters | None
    parameters_vib: ParametersVib | None
    results: VibResults


class PhononSchema(RunSchema):
    """Schema for [quacc.schemas.phonons.summarize_phonopy][]"""

    force_constant: NDArray


class ParametersThermo(TypedDict):
    temperature: float
    pressure: float
    sigma: int
    spin_multiplicity: int
    vib_freqs: list[float]
    vib_energies: list[float]
    n_imag: int


class ThermoResults(TypedDict):
    energy: float
    enthalpy: float
    entropy: float
    gibbs_energy: float
    zpe: float


class ThermoSchema(AtomsSchema):
    parameters_thermo: ParametersThermo
    results: ThermoResults


class VibThermoSchema(VibSchema, ThermoSchema):
    """Schema for [quacc.schemas.ase.summarize_vib_and_thermo][]"""
