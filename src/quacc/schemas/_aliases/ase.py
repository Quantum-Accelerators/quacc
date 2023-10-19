"""Aliases for type hinting `quacc.schemas.ase`"""
from typing import Any, TypedDict

from quacc.schemas._aliases.atoms import AtomsSchema

results = dict[str, Any]  # from atoms.calc.results
parameters = dict[str, Any]  # from atoms.calc.parameters


class RunSchema(AtomsSchema):
    """Schema for `quacc.schemas.ase.summarize_run`"""

    input_atoms: AtomsSchema | None
    nid: str
    dir_name: str
    parameters: parameters
    results: results


class OptSchema(RunSchema):
    """Schema for `quacc.schemas.ase.summarize_opt_run`"""

    fmax: float | None
    parameters_opt: dict[str, Any]  # from Optimizer.todict()
    converged: bool
    nsteps: int
    trajectory: list[AtomsSchema]
    trajectory_results: list[results]


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
    """Schema for `quacc.schemas.ase.summarize_vib_run`"""

    nid: str
    dir_name: str
    parameters: parameters | None
    parameters_vib: ParametersVib | None
    results: VibResults


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
    nid: str
    dir_name: str
    parameters_thermo: ParametersThermo
    results: ThermoResults


class VibThermoSchema(VibSchema, ThermoSchema):
    """Schema for `quacc.schemas.ase.summarize_vib_and_thermo`"""
