"""Common elastic constants workflows."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase import units
from ase.stress import voigt_6_to_full_3x3_stress
from emmet.core.elasticity import ElasticityDoc
from emmet.core.mpid import MPID
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.io.ase import AseAtomsAdaptor

from quacc import job
from quacc.atoms.deformation import make_deformations_from_bulk

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms
    from pymatgen.analysis.elasticity.strain import DeformedStructureSet

    from quacc import Job
    from quacc.types import ElasticSchema, OptSchema, RunSchema


@job
def deformations_to_elastic_tensor(
    undeformed_result: OptSchema | RunSchema,
    deformed_structure_set: DeformedStructureSet,
    results: list[dict],
) -> ElasticityDoc:
    structure = AseAtomsAdaptor.get_structure(undeformed_result["atoms"])  # type: ignore
    return ElasticityDoc.from_deformations_and_stresses(
        structure,
        material_id=MPID("quacc-00"),
        deformations=deformed_structure_set.deformations,
        equilibrium_stress=Stress(
            (
                voigt_6_to_full_3x3_stress(undeformed_result["results"]["stress"])
                if len(undeformed_result["results"]["stress"]) == 6
                else undeformed_result["results"]["stress"]
            )
            / units.GPa
        ),
        stresses=[
            Stress(
                (
                    voigt_6_to_full_3x3_stress(relax_result["results"]["stress"])
                    if len(undeformed_result["results"]["stress"]) == 6
                    else undeformed_result["results"]["stress"]
                )
                / units.GPa
            )
            for relax_result in results
        ],
    )


def bulk_to_deformations_subflow(
    atoms: Atoms,
    relax_job: Job,
    static_job: Job,
    run_static: bool = False,
    pre_relax: bool = True,
    deform_kwargs: dict[str, Any] | None = None,
) -> ElasticSchema:
    """
    Workflow consisting of:

    1. Deformed structures generation

    2. Deformed structures relaxations

    3. Deformed structures statics (optional)

    Parameters
    ----------
    atoms
        Atoms object
    relax_job
        The relaxation function.
    static_job
        The static function
    pre_relax
        Whether to pre-relax the input atoms as is common
    static_job
        The static function.
    deform_kwargs
        Additional keyword arguments to pass to
        [quacc.atoms.deformation.make_deformations_from_bulk][]

    Returns
    -------
    list[dict]
        List of schemas.
    """
    deform_kwargs = deform_kwargs or {}

    if pre_relax:
        undeformed_result = relax_job(atoms, relax_cell=True)
    else:
        undeformed_result = static_job(atoms)

    deformed_structure_set = make_deformations_from_bulk(
        undeformed_result["atoms"], **deform_kwargs
    )

    results = []
    for deformed in deformed_structure_set:
        result = relax_job(deformed.to_ase_atoms())

        if run_static:
            result = static_job(result["atoms"])

        results.append(result)

    elasticity_doc = deformations_to_elastic_tensor(
        undeformed_result=undeformed_result,
        deformed_structure_set=deformed_structure_set,
        results=results,
    )

    return {
        "deformed_structure_set": deformed_structure_set,
        "deformed_results": results,
        "undeformed_result": undeformed_result,
        "elasticity_doc": elasticity_doc,
    }
