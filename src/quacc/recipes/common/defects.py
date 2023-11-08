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


def common_bulk_to_defects_flow(
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
) -> list:
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
    list
        List of dictionary of results
    """
    defect_relax_kwargs = defect_relax_kwargs or {}
    defect_static_kwargs = defect_static_kwargs or {}
    make_defects_kwargs = make_defects_kwargs or {}

    @subflow
    def _relax_job_distributed(atoms: Atoms) -> list:
        defects = make_defects_from_bulk(
            atoms,
            defect_gen=defect_gen,
            defect_charge=defect_charge,
            **make_defects_kwargs,
        )
        return [relax_job(defect, **defect_relax_kwargs) for defect in defects]

    @subflow
    def _relax_and_static_job_distributed(atoms: Atoms) -> list:
        defects = make_defects_from_bulk(
            atoms,
            defect_gen=defect_gen,
            defect_charge=defect_charge,
            **make_defects_kwargs,
        )
        return [
            static_job(
                relax_job(defect, **defect_relax_kwargs)["atoms"],
                **defect_static_kwargs,
            )
            for defect in defects
        ]

    return (
        _relax_and_static_job_distributed(atoms)
        if static_job
        else _relax_job_distributed(atoms)
    )
