"""Common defect workflows."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import subflow
from quacc.atoms.defects import make_defects_from_bulk

has_pmg_defects = bool(find_spec("pymatgen.analysis.defects"))
has_shakenbreak = bool(find_spec("shakenbreak"))


if TYPE_CHECKING:
    from typing import Any, TypeVar

    from ase.atoms import Atoms

    from quacc import Job
    from quacc.types import RunSchema

    ResultSchema = TypeVar("ResultSchema", bound=RunSchema)


@subflow
@requires(
    has_pmg_defects,
    "Missing pymatgen-analysis-defects. Please run pip install quacc[defects]",
)
@requires(has_shakenbreak, "Missing shakenbreak. Please run pip install quacc[defects]")
def bulk_to_defects_subflow(
    atoms: Atoms,
    relax_job: Job[ResultSchema],
    static_job: Job[ResultSchema] | None = None,
    make_defects_kwargs: dict[str, Any] | None = None,
) -> dict[str, list[ResultSchema]]:
    """
    Workflow consisting of:

    1. Defect generation

    2. Defect relaxations

    3. Defect statics (optional)

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    relax_job
        Relaxation job. Its result schema is preserved in the returned lists.
    static_job
        Optional static job. Its result schema is preserved in the returned lists.
    make_defects_kwargs
        Keyword arguments for [quacc.atoms.defects.make_defects_from_bulk][]

    Returns
    -------
    dict[str, list[ResultSchema]]
        Results grouped under "relax" and "static". Each list contains the result
        schema returned by the corresponding job, such as
        [RunSchema][quacc.types.RunSchema] or [VaspSchema][quacc.types.VaspSchema].
        The "static" list is empty when no static job is supplied.
    """
    make_defects_kwargs = make_defects_kwargs or {}
    defects = make_defects_from_bulk(atoms, **make_defects_kwargs)

    results: dict[str, list[ResultSchema]] = {"relax": [], "static": []}
    for defect in defects:
        result = relax_job(defect)
        results["relax"].append(result)

        if static_job is not None:
            result = static_job(result["atoms"])
            results["static"].append(result)

    return results
