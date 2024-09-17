"""Defect recipes for EMT."""

from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import flow
from quacc.recipes.common.defects import bulk_to_defects_subflow
from quacc.recipes.emt.core import relax_job, static_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

has_pmg_defects = bool(find_spec("pymatgen.analysis.defects"))
has_shakenbreak = bool(find_spec("shakenbreak"))

if has_pmg_defects:
    from pymatgen.analysis.defects.generators import VacancyGenerator


if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import OptSchema, RunSchema

    if has_pmg_defects:
        from pymatgen.analysis.defects.generators import (
            AntiSiteGenerator,
            ChargeInterstitialGenerator,
            InterstitialGenerator,
            SubstitutionGenerator,
            VoronoiInterstitialGenerator,
        )


@flow
@requires(
    has_pmg_defects,
    "Missing pymatgen-analysis-defects. Please run pip install quacc[defects]",
)
@requires(has_shakenbreak, "Missing shakenbreak. Please run pip install quacc[defects]")
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
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> list[RunSchema | OptSchema]:
    """
    Workflow consisting of:

    1. Defect generation

    2. Defect relaxations
        - name: "relax_job"
        - job: [quacc.recipes.emt.core.relax_job][]

    3. Optional defect statics
        - name: "static_job"
        - job: [quacc.recipes.emt.core.static_job][]

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
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictinoary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    list[RunSchema | OptSchema]
        List of dictionary of results from [quacc.schemas.ase.Summarize.run][]
        or [quacc.schemas.ase.Summarize.opt][].
        See the return type-hint for the data structure.
    """
    relax_job_, static_job_ = customize_funcs(
        ["relax_job", "static_job"],
        [relax_job, static_job],
        param_swaps=job_params,
        decorators=job_decorators,
    )
    make_defects_kwargs = recursive_dict_merge(
        make_defects_kwargs, {"defect_gen": defect_gen, "defect_charge": defect_charge}
    )

    return bulk_to_defects_subflow(
        atoms,
        relax_job_,
        static_job=static_job_ if run_static else None,
        make_defects_kwargs=make_defects_kwargs,
    )
