"""Aliases for type hinting `quacc.schemas.ase`"""
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any, TypedDict

    from quacc.schemas._aliases.atoms import AtomsSchema

    results = dict[str, Any]  # from atoms.calc.results

    class RunSchema(TypedDict):
        """Schema for `summarize_run`"""

        input_atoms: AtomsSchema | None
        nid: str
        dir_name: str
        parameters: dict[str, Any]  # from atoms.calc.parameters
        results: results

    class OptSchema(RunSchema):
        """Schema for `summarize_opt_run`"""

        fmax: float | None
        parameters_opt: dict[str, Any]  # from Optimizer.todict()
        converged: bool
        nsteps: int
        trajectory: list[AtomsSchema]
        trajectory_results: list[results]

    class VibInputs(TypedDict):
        dir_name: str
        parameters: dict[str, Any]
        paramters_vib: dict[str, Any]
        nid: str

    class VibDataInputs(TypedDict):
        dir_name: str
        delta: float
        direction: str
        method: str
        ndof: int
        nfree: int
        nid: str

    class VibResults(TypedDict):
        imag_vib_freqs: list[float]
        n_imag: int
        vib_energies: list[float]
        vib_freqs: list[float]
        vib_energies_raw: list[float]
        vib_freqs_raw: list[float]

    class VibSchema(TypedDict):
        """Schema for `summarize_vib_run`"""

        parameters_vib: dict[str, Any]
        results: VibResults

    class ParamatersThermo(TypedDict):
        temperature: float
        pressure: float
        sigma: float
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

    class ThermoSchema(BaseRunnerSchema):
        parameters_thermo: ParamatersThermo
        results: ThermoResults

    class FreqSchema(VibSchema):
        thermo: ThermoSchema
