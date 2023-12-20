"""Defect recipes for EMT."""
from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.analysis.defects.generators import VacancyGenerator

from quacc import flow
from quacc.recipes.common.defects import bulk_to_defects_subflow
from quacc.recipes.emt.core import relax_job, static_job
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from pymatgen.analysis.defects.generators import (
        AntiSiteGenerator,
        ChargeInterstitialGenerator,
        InterstitialGenerator,
        SubstitutionGenerator,
        VoronoiInterstitialGenerator,
    )

    from quacc import Job
    from quacc.schemas._aliases.ase import OptSchema, RunSchema


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
    custom_relax_job: Job | None = None,
    custom_static_job: Job | None = None,
    run_static: bool = True,
    make_defects_kwargs: dict[str, Any] | None = None,
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
    custom_relax_job
        Relaxation job, which defaults to [quacc.recipes.emt.core.relax_job][].
    custom_static_job
        Static job, which defaults to [quacc.recipes.emt.core.static_job][].
    make_defects_kwargs
        Keyword arguments to pass to
        [quacc.atoms.defects.make_defects_from_bulk][]

    Returns
    -------
    list[RunSchema | OptSchema]
        List of dictionary of results from [quacc.schemas.ase.summarize_run][]
        or [quacc.schemas.ase.summarize_opt_run][]
    """
    make_defects_kwargs = merge_dicts(
        make_defects_kwargs, {"defect_gen": defect_gen, "defect_charge": defect_charge}
    )

    return bulk_to_defects_subflow(
        atoms,
        relax_job if custom_relax_job is None else custom_relax_job,
        static_job=(static_job if custom_static_job is None else custom_static_job)
        if run_static
        else None,
        make_defects_kwargs=make_defects_kwargs,
    )
