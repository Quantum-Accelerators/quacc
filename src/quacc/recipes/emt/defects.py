"""Defect recipes for EMT"""
from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.analysis.defects.generators import VacancyGenerator

from quacc import flow, job, subflow
from quacc.recipes.emt.core import relax_job as _relax_job
from quacc.recipes.emt.core import static_job as _static_job
from quacc.schemas import fetch_atoms
from quacc.utils.defects import make_defects_from_bulk

if TYPE_CHECKING:
    from ase import Atoms
    from pymatgen.analysis.defects.generators import (
        AntiSiteGenerator,
        ChargeInterstitialGenerator,
        InterstitialGenerator,
        SubstitutionGenerator,
        VoronoiInterstitialGenerator,
    )

    from quacc.schemas.ase import OptSchema, RunSchema
    from quacc.utils.wflows import Job


@flow
def bulk_to_defects_flow(
    atoms: Atoms | dict,
    defect_gen: (
        AntiSiteGenerator
        | ChargeInterstitialGenerator
        | InterstitialGenerator
        | SubstitutionGenerator
        | VacancyGenerator
        | VoronoiInterstitialGenerator
    ) = VacancyGenerator,
    defect_charge: int = 0,
    make_defects_kwargs: dict | None = None,
    defect_relax: Job | None = _relax_job,
    defect_static: Job | None = _static_job,
    defect_relax_kwargs: dict | None = None,
    defect_static_kwargs: dict | None = None,
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
        Keyword arguments to pass to the make_defects_from_bulk
    defect_relax
        Default Job to use for the relaxation of the defect structures.
    defect_static
        Default Job to use for the static calculation of the defect structures.
    defect_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    defect_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    list[dict]
        List of dictionary of results from `quacc.schemas.ase.summarize_run` or
        `quacc.schemas.ase.summarize_opt_run`
    """
    defect_relax_kwargs = defect_relax_kwargs or {}
    defect_static_kwargs = defect_static_kwargs or {}
    make_defects_kwargs = make_defects_kwargs or {}

    if "relax_cell" not in defect_relax_kwargs:
        defect_relax_kwargs["relax_cell"] = False

    @job
    def _make_defects(atoms):
        atoms = fetch_atoms(atoms)
        return make_defects_from_bulk(
            atoms,
            defect_gen=defect_gen,
            defect_charge=defect_charge,
            **make_defects_kwargs,
        )

    @subflow
    def _relax_distributed(defects):
        return [defect_relax(defect, **defect_relax_kwargs) for defect in defects]

    @subflow
    def _relax_and_static_distributed(defects):
        return [
            defect_static(
                defect_relax(defect, **defect_relax_kwargs),
                **defect_static_kwargs,
            )
            for defect in defects
        ]

    defects = _make_defects(atoms)

    if defect_static is None:
        return _relax_distributed(defects)

    return _relax_and_static_distributed(defects)
