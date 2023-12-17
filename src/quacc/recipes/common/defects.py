"""Common defect workflows"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import subflow
from quacc.atoms.defects import make_defects_from_bulk

if TYPE_CHECKING:
    from typing import Callable

    from ase.atoms import Atoms

    from quacc import Job


@subflow
def bulk_to_defects_subflow(
    atoms: Atoms,
    relax_job: Job,
    static_job: Job | None = None,
    make_defects_fn: Callable = make_defects_from_bulk,
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
    make_defects_fn
        The function for generating defects.

    Returns
    -------
    list[dict]
        List of dictionary of results
    """

    defects = make_defects_fn(atoms)

    results = []
    for defect in defects:
        result = relax_job(defect)

        if static_job:
            result = static_job(result["atoms"])

        results.append(result)

    return results
