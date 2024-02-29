"""Common elastic constants workflows."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import subflow
from quacc.atoms.deformation import make_deformations_from_bulk

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc import Job


@subflow
def bulk_to_deformations_subflow(
    atoms: Atoms,
    relax_job: Job,
    static_job: Job | None = None,
    deform_kwargs: dict[str, Any] | None = None,
) -> list[dict]:
    """
    Workflow consisting of:

    1. Deformed structures generation

    2. Deformed structures relaxations

    3. Deformed structures statics (optional)

    Parameters
    ----------
    atoms
        Atoms object
    relax_job
        The relaxation function.
    static_job
        The static function.
    deform_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.deformation.make_deformations_from_bulk][]

    Returns
    -------
    list[dict]
        List of schemas.
    """
    deform_kwargs = deform_kwargs or {}

    deformations = make_deformations_from_bulk(atoms, **deform_kwargs)

    results = []
    for deformed in deformations:
        result = relax_job(deformed)

        if static_job is not None:
            result = static_job(result["atoms"])

        results.append(result)

    return results
