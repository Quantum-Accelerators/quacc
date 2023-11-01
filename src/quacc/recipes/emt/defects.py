"""Defect recipes for EMT."""
from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.analysis.defects.generators import VacancyGenerator

from quacc import flow, job, subflow
from quacc.atoms.defects import make_defects_from_bulk
from quacc.recipes.emt.core import relax_job, static_job

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms
    from pymatgen.analysis.defects.generators import (
        AntiSiteGenerator,
        ChargeInterstitialGenerator,
        InterstitialGenerator,
        SubstitutionGenerator,
        VoronoiInterstitialGenerator,
    )

    from quacc.schemas.ase import OptSchema, RunSchema


@flow
def bulk_to_defects_flow(
    atoms: Atoms,
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
    run_static: bool = True,
    defect_relax_kwargs: dict[str, Any] | None = None,
    defect_static_kwargs: dict[str, Any] | None = None,
) -> list[RunSchema | OptSchema]:
    """
    Workflow consisting of:

    1. Defect generation

    2. Defect relaxations

    3. Defect statics (optional)

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    defect_gen
        Defect generator
    defect_charge
        Charge state of the defect
    make_defects_kwargs
        Keyword arguments to pass to
        [quacc.atoms.defects.make_defects_from_bulk][]
    run_static
        Whether to run the static calculation.
    defect_relax_kwargs
        Additional keyword arguments to pass to [quacc.recipes.emt.core.relax_job][].
    defect_static_kwargs
        Additional keyword arguments to pass to [quacc.recipes.emt.core.static_job][].

    Returns
    -------
    list[RunSchema | OptSchema]
        List of dictionary of results from [quacc.schemas.ase.summarize_run][]
        or [quacc.schemas.ase.summarize_opt_run][]
    """
    defect_relax_kwargs = defect_relax_kwargs or {}
    defect_static_kwargs = defect_static_kwargs or {}
    make_defects_kwargs = make_defects_kwargs or {}

    if "relax_cell" not in defect_relax_kwargs:
        defect_relax_kwargs["relax_cell"] = False

    @job
    def _make_defects(atoms):
        return make_defects_from_bulk(
            atoms,
            defect_gen=defect_gen,
            defect_charge=defect_charge,
            **make_defects_kwargs,
        )

    @subflow
    def _relax_distributed(defects):
        return [relax_job(defect, **defect_relax_kwargs) for defect in defects]

    @subflow
    def _relax_and_static_distributed(defects):
        return [
            static_job(
                relax_job(defect, **defect_relax_kwargs)["atoms"],
                **defect_static_kwargs,
            )
            for defect in defects
        ]

    defects = _make_defects(atoms)

    return (
        _relax_and_static_distributed(defects)
        if run_static
        else _relax_distributed(defects)
    )
