"""Elastic constants recipes for EMT."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.common.elastic import elastic_tensor_flow as elastic_tensor_flow_
from quacc.recipes.emt.core import relax_job, static_job
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import ElasticSchema


@flow
def elastic_tensor_flow(
    atoms: Atoms,
    pre_relax: bool = True,
    run_static: bool = False,
    deform_kwargs: dict[str, Any] | None = None,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> ElasticSchema:
    """
    Workflow consisting of:

    1. Bulk structure relaxation (if pre_relax is True)
        - name: "relax_job"
        - job: [quacc.recipes.emt.core.relax_job][]

    2. Bulk structure static calculation (if run_static is True)
        - name: "static_job"
        - job: [quacc.recipes.emt.core.static_job][]

    3. Deformed structures generation

    4. Deformed structures relaxations
        - name: "relax_job"
        - job: [quacc.recipes.emt.core.relax_job][]

    5. Deformed structures statics (if run_static is True)
        - name: "static_job"
        - job: [quacc.recipes.emt.core.static_job][]

    6. Elastic tensor calculation

    Parameters
    ----------
    atoms
        Atoms object
    pre_relax
        Whether to run a relaxation on the structure before deformation (true)
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

    return elastic_tensor_flow_(
        atoms=atoms,
        relax_job=relax_job_,
        static_job=static_job_,
        pre_relax=pre_relax,
        run_static=run_static,
        deform_kwargs=deform_kwargs,
    )
