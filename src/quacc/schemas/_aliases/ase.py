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
