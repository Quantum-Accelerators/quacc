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
    from typing import Any

    from ase.atoms import Atoms

    from quacc import Job


@subflow
@requires(
    has_pmg_defects,
    "Missing pymatgen-analysis-defects. Please run pip install quacc[defects]",
)
@requires(has_shakenbreak, "Missing shakenbreak. Please run pip install quacc[defects]")
def bulk_to_defects_subflow(
    atoms: Atoms,
    relax_job: Job,
    static_job: Job | None = None,
    make_defects_kwargs: dict[str, Any] | None = None,
) -> list[dict]:
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
        The relaxation function.
    static_job
        The static function.
    make_defects_kwargs
        Keyword arguments for [quacc.atoms.defects.make_defects_from_bulk][]

    Returns
    -------
    list[dict]
        List of dictionary of results.
    """
    make_defects_kwargs = make_defects_kwargs or {}
    defects = make_defects_from_bulk(atoms, **make_defects_kwargs)

    results = []
    for defect in defects:
        result = relax_job(defect)

        if static_job is not None:
            result = static_job(result["atoms"])

        results.append(result)

    return results
