"""Common defect workflows"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import subflow
from quacc.atoms.defects import make_defects_from_bulk

if TYPE_CHECKING:
    from typing import Callable

    from ase import Atoms


@subflow
def bulk_to_defects_subflow(
    atoms: Atoms,
    relax_fn: Callable,
    static_fn: Callable | None = None,
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
    relax_fn
        The relaxation function.
    static_fn
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
        result = relax_fn(defect)

        if static_fn:
            result = static_fn(result["atoms"])

        results.append(result)

    return results
