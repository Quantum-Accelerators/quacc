"""Elastic constants recipes for MLPs."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.common.elastic import bulk_to_deformations_subflow
from quacc.recipes.mlp.core import relax_job, static_job
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import ElasticSchema


@flow
def bulk_to_deformations_flow(
    atoms: Atoms,
    pre_relax: bool = True,
    run_static: bool = False,
    deform_kwargs: dict[str, Any] | None = None,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> ElasticSchema:
    """
    Workflow consisting of:

    1. Deformed structures generation

    2. Deformed structures relaxations
        - name: "relax_job"
        - job: [quacc.recipes.mlp.core.relax_job][]

    3. Deformed structures statics (optional)
        - name: "static_job"
        - job: [quacc.recipes.mlp.core.static_job][]

    Parameters
    ----------
    atoms
        Atoms object
    pre_relax
        Whether to pre-relax the input atoms as is common
    run_static
        Whether to run static calculations after the relaxations
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
    )  # type: ignore

    if pre_relax:
        if run_static:
            undeformed_relax = relax_job_(atoms, relax_cell=True)
            undeformed_result = static_job_(undeformed_relax["atoms"])
        else:
            undeformed_result = relax_job_(atoms, relax_cell=True)
    else:
        undeformed_result = static_job_(atoms)

    return bulk_to_deformations_subflow(
        undeformed_result=undeformed_result,
        relax_job=relax_job_,
        static_job=static_job_,
        run_static=run_static,
        deform_kwargs=deform_kwargs,
    )
