"""Defect recipes for EMT."""
from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.analysis.defects.generators import VacancyGenerator

from quacc import flow
from quacc.recipes.common.defects import bulk_to_defects_subflow
from quacc.recipes.emt.core import relax_job, static_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms
    from pymatgen.analysis.defects.generators import (
        AntiSiteGenerator,
        ChargeInterstitialGenerator,
        InterstitialGenerator,
        SubstitutionGenerator,
        VoronoiInterstitialGenerator,
    )

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
    run_static: bool = True,
    make_defects_kwargs: dict[str, Any] | None = None,
    decorators: dict[str, Callable | None] | None = None,
    parameters: dict[str, Any] | None = None,
) -> list[RunSchema | OptSchema]:
    """
    Workflow consisting of:

    1. Defect generation

    2. Defect relaxations ("relax_job")

    3. Optional defect statics ("static_job")

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    defect_gen
        Defect generator
    defect_charge
        Charge state of the defect
    run_static
        Whether to run static calculations.
    make_defects_kwargs
        Keyword arguments to pass to
        [quacc.atoms.defects.make_defects_from_bulk][]
    decorators
        Custom decorators to apply to each Job in the Flow.
        Refer to [quacc.wflow_tools.customizers.customize_funcs][] for details.
    parameters
        Custom parameters to pass to each Job in the Flow.
        Refer to [quacc.wflow_tools.customizers.customize_funcs][] for details.

    Returns
    -------
    list[RunSchema | OptSchema]
        List of dictionary of results from [quacc.schemas.ase.summarize_run][]
        or [quacc.schemas.ase.summarize_opt_run][]
    """
    make_defects_kwargs = recursive_dict_merge(
        make_defects_kwargs, {"defect_gen": defect_gen, "defect_charge": defect_charge}
    )
    relax_job_, static_job_ = customize_funcs(
        {"relax_job": relax_job, "static_job": static_job},
        decorators=decorators,
        parameters=parameters,
    )

    return bulk_to_defects_subflow(
        atoms,
        relax_job_,
        static_job=static_job_ if run_static else None,
        make_defects_kwargs=make_defects_kwargs,
    )
