"""Common elastic constants workflows."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase import units
from ase.stress import voigt_6_to_full_3x3_stress
from emmet.core.elasticity import ElasticityDoc
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.io.ase import AseAtomsAdaptor

from quacc import flow, job, subflow
from quacc.atoms.deformation import make_deformations_from_bulk
from quacc.wflow_tools.customizers import customize_jobs

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from ase.atoms import Atoms
    from pymatgen.analysis.elasticity.strain import DeformedStructureSet

    from quacc import Job
    from quacc.types import ElasticSchema, OptSchema, RunSchema


@flow
def elastic_tensor_flow(
    atoms: Atoms,
    relax_job: Job,
    static_job: Job,
    pre_relax: bool = True,
    run_static: bool = False,
    deform_kwargs: dict[str, Any] | None = None,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> ElasticSchema:
    """
    Common workflow for calculating elastic tensors.

    Parameters
    ----------
    atoms
        Atoms object
    relax_job
        The relaxation function.
    static_job
        The static function
    pre_relax
        Whether to run a relaxation on the bulk structure before deformation (true) or run a static
        calculation (false). The bulk cell is relaxed by default. Set
        `job_params={"relax_job": {"relax_cell": False}}` to keep it fixed.
        Deformed structures are always relaxed at fixed cell.
    run_static
        Whether to run static calculations after any relaxations on the undeformed or deformed structures
    deform_kwargs
        Additional keyword arguments to pass to [quacc.atoms.deformation.make_deformations_from_bulk][]
    job_params
        Custom parameters for `relax_job` and `static_job`, keyed by job name.
        The key `"all"` applies parameters to both jobs before job-specific
        parameters. The relaxation job defaults to `relax_cell=True`; deformed
        structures always override this with `relax_cell=False`.
    job_decorators
        Custom decorators for `relax_job` and `static_job`, keyed by job name.
        The key `"all"` applies a decorator to both jobs. A value of None strips
        the existing decorator.

    Returns
    -------
    ElasticSchema
        See the return type-hint for the data structure.
    """
    relax_job, static_job = customize_jobs(
        {"relax_job": relax_job, "static_job": static_job},
        param_defaults={"relax_job": {"relax_cell": True}},
        param_swaps=job_params,
        decorators=job_decorators,
    ).values()

    if pre_relax:
        undeformed_result = relax_job(atoms)
        if run_static:
            undeformed_result = static_job(undeformed_result["atoms"])
    else:
        undeformed_result = static_job(atoms)

    return _elastic_tensor_subflow(
        undeformed_result=undeformed_result,
        relax_job=relax_job,
        static_job=static_job if run_static else None,
        deform_kwargs=deform_kwargs,
    )


@subflow
def _elastic_tensor_subflow(
    undeformed_result: RunSchema | OptSchema,
    relax_job: Job,
    static_job: Job | None = None,
    deform_kwargs: dict[str, Any] | None = None,
) -> ElasticSchema:
    """
    Workflow consisting of:

    1. Deformed structures generation

    2. Deformed structures relaxations

    3. Deformed structures statics (optional)

    4. Elastic tensor calculation

    Parameters
    ----------
    undeformed_result
        Result of a static or optimization calculation
    relax_job
        The relaxation function.
    static_job
        The static function
    deform_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.deformation.make_deformations_from_bulk][]

    Returns
    -------
    ElasticSchema
        Elastic properties with deformed results grouped under "relax" and "static".
        The "static" list is empty when no static job is supplied.
    """
    deform_kwargs = deform_kwargs or {}

    deformed_structure_set = make_deformations_from_bulk(
        undeformed_result["atoms"], **deform_kwargs
    )

    results: dict[str, list[RunSchema | OptSchema]] = {"relax": [], "static": []}
    for deformed in deformed_structure_set:
        result = relax_job(deformed.to_ase_atoms(), relax_cell=False)
        results["relax"].append(result)

        if static_job is not None:
            result = static_job(result["atoms"])
            results["static"].append(result)

    elasticity_doc = _deformations_to_elastic_tensor(
        undeformed_result=undeformed_result,
        deformed_structure_set=deformed_structure_set,
        results=results["relax"],
    )

    return {
        "deformed_structure_set": deformed_structure_set,
        "deformed_results": results,
        "undeformed_result": undeformed_result,
        "elasticity_doc": elasticity_doc,
    }


@job
def _deformations_to_elastic_tensor(
    undeformed_result: RunSchema | OptSchema,
    deformed_structure_set: DeformedStructureSet,
    results: list[RunSchema | OptSchema],
) -> ElasticityDoc:
    """
    Function to fit a DeformedStructureSet and result documents to an elastic tensor

    Parameters
    ----------
    undeformed_result
        Single point or relaxation of the undeformed result (for the relaxed stress, if it's not quite 0)
    deformed_structure_set
        The pymatgen DeformedStructureSet with information on the strains of each structure
    results
        A list of results (one per deformed structure) corresponding to single points or relaxations
        on the deformed structures.

    Returns
    -------
    ElasticityDoc
        An emmet (i.e. MaterialsProject) ElasticityDoc with information about the final elastic tensor fit
    """
    structure = AseAtomsAdaptor.get_structure(undeformed_result["atoms"])
    equilibrium_stress = Stress(
        (
            voigt_6_to_full_3x3_stress(undeformed_result["results"]["stress"])
            if len(undeformed_result["results"]["stress"]) == 6
            else undeformed_result["results"]["stress"]
        )
        / units.GPa
    )
    stresses = [
        Stress(
            (
                voigt_6_to_full_3x3_stress(relax_result["results"]["stress"])
                if len(relax_result["results"]["stress"]) == 6
                else relax_result["results"]["stress"]
            )
            / units.GPa
        )
        for relax_result in results
    ]
    return ElasticityDoc.from_deformations_and_stresses(
        structure,
        deformations=deformed_structure_set.deformations,
        equilibrium_stress=equilibrium_stress,
        stresses=stresses,
    )
