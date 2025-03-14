"""Elastic constants recipes for MLPs."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.common.elastic import bulk_to_elastic_tensor_subflow
from quacc.recipes.mlp.core import relax_job, static_job
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import ElasticSchema


@flow
def bulk_to_elastic_tensor_flow(
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
        - job: [quacc.recipes.emt.core.relax_job][]

    3. Deformed structures statics (optional)
        - name: "static_job"
        - job: [quacc.recipes.emt.core.static_job][]

    4. Elastic tensor calculation


    Parameters
    ----------
    atoms
        Atoms object
    pre_relax
        Whether to run a relaxation on the bulk structure before deformation (true) or run a static
        calculation (false)
    run_static
        Whether to run static calculations after any relaxations on the undeformed or deformed structures
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
    ElasticSchema
        See the return type-hint for the data structure.
    """
    relax_job_, static_job_ = customize_funcs(
        ["relax_job", "static_job"],
        [relax_job, static_job],
        param_swaps=job_params,
        decorators=job_decorators,
    )  # type: ignore

    if pre_relax:
        undeformed_result = relax_job_(atoms, relax_cell=True)
        if run_static:
            undeformed_result = static_job_(undeformed_result["atoms"])
    else:
        undeformed_result = static_job_(atoms)

    return bulk_to_elastic_tensor_subflow(
        undeformed_result=undeformed_result,
        relax_job=relax_job_,
        static_job=static_job_ if run_static else None,
        deform_kwargs=deform_kwargs,
    )
