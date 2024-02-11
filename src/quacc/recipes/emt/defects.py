"""Defect recipes for EMT."""

from __future__ import annotations

from typing import TYPE_CHECKING

from monty.dev import requires

from quacc import flow
from quacc.recipes.common.defects import bulk_to_defects_subflow
from quacc.recipes.emt.core import relax_job, static_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

try:
    import shakenbreak  # noqa: F401
    from pymatgen.analysis.defects.generators import VacancyGenerator

    has_deps = True
except ImportError:
    has_deps = False


if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import OptSchema, RunSchema

    if has_deps:
        from pymatgen.analysis.defects.generators import (
            AntiSiteGenerator,
            ChargeInterstitialGenerator,
            InterstitialGenerator,
            SubstitutionGenerator,
            VoronoiInterstitialGenerator,
        )


@flow
@requires(
    has_deps, "Missing defect dependencies. Please run pip install quacc[defects]"
)
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
        List of dictionary of results from [quacc.schemas.ase.summarize_run][]
        or [quacc.schemas.ase.summarize_opt_run][].
        See the return type-hint for the data structure.
    """
    make_defects_kwargs = recursive_dict_merge(
        make_defects_kwargs, {"defect_gen": defect_gen, "defect_charge": defect_charge}
    )
    relax_job_, static_job_ = customize_funcs(
        ["relax_job", "static_job"],
        [relax_job, static_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    return bulk_to_defects_subflow(
        atoms,
        relax_job_,
        static_job=static_job_ if run_static else None,
        make_defects_kwargs=make_defects_kwargs,
    )
