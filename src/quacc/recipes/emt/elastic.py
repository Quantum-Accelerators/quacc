"""Elastic constants recipes for EMT."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.common.elastic import bulk_to_deformations_subflow
from quacc.recipes.emt.core import relax_job, static_job
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable

    from ase.atoms import Atoms

    from quacc.types import OptSchema, RunSchema


@flow
def bulk_to_deformations_flow(
    atoms: Atoms,
    run_static: bool = True,
    deform_kwargs: dict[str, Any] | None = None,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> list[RunSchema | OptSchema]:
    """
    Workflow consisting of:

    1. Deformed structures generation

    2. Deformed structures relaxations
        - name: "relax_job"
        - job: [quacc.recipes.emt.core.relax_job][]

    3. Deformed structures statics (optional)
        - name: "static_job"
        - job: [quacc.recipes.emt.core.static_job][]

    Parameters
    ----------
    atoms
        Atoms object
    run_static
        Whether to run static calculations.
    deform_kwargs
        Additional keyword arguments to pass to [quacc.atoms.deformation.make_deformations_from_bulk][]
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    list[RunSchema | OptSchema]
        [RunSchema][quacc.schemas.ase.Summarize.run] or
        [OptSchema][quacc.schemas.ase.Summarize.opt] for each deformation.
        See the return type-hint for the data structure.
    """
    relax_job_, static_job_ = customize_funcs(
        ["relax_job", "static_job"],
        [relax_job, static_job],
        param_swaps=job_params,
        decorators=job_decorators,
    )

    return bulk_to_deformations_subflow(
        atoms,
        relax_job_,
        static_job=static_job_ if run_static else None,
        deform_kwargs=deform_kwargs,
    )
