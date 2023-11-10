"""Common defect workflows"""
from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.analysis.defects.generators import VacancyGenerator

from quacc import subflow
from quacc.atoms.defects import make_defects_from_bulk

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase import Atoms
    from pymatgen.analysis.defects.generators import (
        AntiSiteGenerator,
        ChargeInterstitialGenerator,
        InterstitialGenerator,
        SubstitutionGenerator,
        VoronoiInterstitialGenerator,
    )


@subflow
def bulk_to_defects_subflow(
    atoms: Atoms,
    relax_job: Callable,
    static_job: Callable | None,
    defect_gen: (
        AntiSiteGenerator
        | ChargeInterstitialGenerator
        | InterstitialGenerator
        | SubstitutionGenerator
        | VacancyGenerator
        | VoronoiInterstitialGenerator
    ) = VacancyGenerator,
    defect_charge: int = 0,
    make_defects_kwargs: dict[str, Any] | None = None,
    defect_relax_kwargs: dict[str, Any] | None = None,
    defect_static_kwargs: dict[str, Any] | None = None,
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
    defect_gen
        Defect generator
    defect_charge
        Charge state of the defect
    make_defects_kwargs
        Keyword arguments to pass to
        [quacc.atoms.defects.make_defects_from_bulk][]
    defect_relax_kwargs
        Additional keyword arguments to pass to [quacc.recipes.emt.core.relax_job][].
    defect_static_kwargs
        Additional keyword arguments to pass to [quacc.recipes.emt.core.static_job][].

    Returns
    -------
    list[dict]
        List of dictionary of results
    """
    defect_relax_kwargs = defect_relax_kwargs or {}
    defect_static_kwargs = defect_static_kwargs or {}
    make_defects_kwargs = make_defects_kwargs or {}

    defects = make_defects_from_bulk(
        atoms,
        defect_gen=defect_gen,
        defect_charge=defect_charge,
        **make_defects_kwargs,
    )

    results = []
    for defect in defects:
        result = relax_job(defect, **defect_relax_kwargs)

        if static_job:
            result = static_job(defect, **defect_static_kwargs)

        results.append(result)

    return results
