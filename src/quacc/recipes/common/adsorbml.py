from __future__ import annotations

import copy
import logging
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
from ase.atoms import Atoms
from ase.build import molecule
from fairchem.data.oc.core.adsorbate import Adsorbate
from fairchem.data.oc.core.bulk import Bulk
from fairchem.data.oc.core.multi_adsorbate_slab_config import (
    MultipleAdsorbateSlabConfig,
)
from fairchem.data.oc.core.slab import Slab, tile_and_tag_atoms
from fairchem.data.oc.utils import DetectTrajAnomaly

from quacc import Job, flow, job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs, strip_decorator

if TYPE_CHECKING:
    from quacc.types import OptSchema, RunSchema

logger = logging.getLogger(__name__)


@job
def ocp_surface_generator(bulk_atoms, max_miller: int = 1):
    # This gets ocdata.core.slab.Slab objects which already know which atoms are the surface atoms!
    return Slab.from_bulk_get_all_slabs(Bulk(bulk_atoms), max_miller)


class CustomSlab(Slab):
    # Custom slab object to ignore the composition of the bulk in the initialization
    # and make sure that the user did the right things to have the surface tagged for
    # adsorbate placement

    def __init__(
        self,
        slab_atoms: Atoms,
        bulk: Atoms | None = None,
        millers: tuple | None = None,
        shift: float | None = None,
        top: bool | None = None,
        min_ab: float = 0.8,
    ):
        self.bulk = bulk
        self.slab_atoms = slab_atoms
        self.millers = millers
        self.shift = shift
        self.top = top

        assert np.linalg.norm(self.slab_atoms.cell[0]) >= min_ab, "Slab not tiled"
        assert np.linalg.norm(self.slab_atoms.cell[1]) >= min_ab, "Slab not tiled"
        assert self.has_surface_tagged(), "Slab not tagged"
        assert len(self.slab_atoms.constraints) > 0, "Sub-surface atoms not constrained"


# Not clear this actually needs to a job, maybe just a function if pretty fast?
@job
def ocp_adslab_generator(
    slab: Slab | Atoms,
    adsorbates_kwargs: list[dict[Any, Any]] | None = None,
    multiple_adsorbate_slab_config_kwargs: dict[Any, Any] | None = None,
) -> list[Atoms]:
    adsorbates = [
        Adsorbate(**adsorbate_kwargs) for adsorbate_kwargs in adsorbates_kwargs
    ]

    if isinstance(slab, Atoms):
        try:
            slab = CustomSlab(slab_atoms=slab)
        except AssertionError:
            slab = CustomSlab(slab_atoms=tile_and_tag_atoms(slab))
            logger.warning(
                "The slab was not tagged and/or tiled. "
                "We did the best we could, but you should be careful and check the results!"
            )

    if multiple_adsorbate_slab_config_kwargs is None:
        multiple_adsorbate_slab_config_kwargs = {}

    adslabs = MultipleAdsorbateSlabConfig(
        copy.deepcopy(slab), adsorbates, **multiple_adsorbate_slab_config_kwargs
    )

    return adslabs.atoms_list


@flow
def find_adslabs_each_slab(
    slabs: list[Slab],
    adsorbates_kwargs: dict[Any, Any],
    multiple_adsorbate_slab_config_kwargs: dict[Any, Any] | None = None,
):
    return [
        {
            "slab": slab,
            "adslabs": ocp_adslab_generator(
                slab, adsorbates_kwargs, multiple_adsorbate_slab_config_kwargs
            ),
        }
        for slab in slabs
    ]


def detect_anomaly(
    initial_atoms: Atoms, final_atoms: Atoms
) -> list[
    Literal[
        "adsorbate_dissociated",
        "adsorbate_desorbed",
        "surface_changed",
        "adsorbate_intercalated",
    ]
]:
    atom_tags = initial_atoms.get_tags()

    detector = DetectTrajAnomaly(initial_atoms, final_atoms, atom_tags)
    anomalies = []
    if detector.is_adsorbate_dissociated():
        anomalies.append("adsorbate_dissociated")
    if detector.is_adsorbate_desorbed():
        anomalies.append("adsorbate_desorbed")
    if detector.has_surface_changed():
        anomalies.append("surface_changed")
    if detector.is_adsorbate_intercalated():
        anomalies.append("adsorbate_intercalated")
    return anomalies


@job
def filter_sort_select_adslabs(
    adslab_results: list[RunSchema], adslab_anomalies_list: list[list[str]]
):
    for adslab_result, adslab_anomalies in zip(
        adslab_results, adslab_anomalies_list, strict=True
    ):
        adslab_result["results"]["adslab_anomalies"] = adslab_anomalies

    adslabs_no_anomalies = [
        adslab_result
        for adslab_result in adslab_results
        if len(adslab_result["results"]["adslab_anomalies"]) == 0
    ]

    return sorted(adslabs_no_anomalies, key=lambda x: x["results"]["energy"])


@flow
def adsorb_ml_pipeline(
    slab: Slab,
    adsorbates_kwargs: dict[Any, Any],
    multiple_adsorbate_slab_config_kwargs: dict[Any, Any],
    ml_slab_adslab_relax_job: Job,
    slab_validate_job: Job,
    adslab_validate_job: Job,
    gas_validate_job: Job,
    num_to_validate_with_DFT: int = 0,
    ml_gas_references: dict[str, RunSchema] | None = None,
    reference_ml_energies_to_gas_phase: bool = False,
):
    unrelaxed_adslab_configurations = ocp_adslab_generator(
        slab, adsorbates_kwargs, multiple_adsorbate_slab_config_kwargs
    )

    ml_relaxed_configurations = [
        ml_slab_adslab_relax_job(adslab_configuration)
        for adslab_configuration in unrelaxed_adslab_configurations
    ]

    ml_relaxed_slab_result = ml_slab_adslab_relax_job(slab.atoms)

    if reference_ml_energies_to_gas_phase:
        if ml_gas_references is None:
            ml_gas_references = {
                "N2": ml_slab_adslab_relax_job(
                    molecule_pbc("N2", vacuum=10), relax_cell=False
                ),
                "CO": ml_slab_adslab_relax_job(
                    molecule_pbc("CO", vacuum=10), relax_cell=False
                ),
                "H2": ml_slab_adslab_relax_job(
                    molecule_pbc("H2", vacuum=10), relax_cell=False
                ),
                "H2O": ml_slab_adslab_relax_job(
                    molecule_pbc("H2O", vacuum=10), relax_cell=False
                ),
            }
        ml_relaxed_configurations = reference_adslab_energies(
            ml_relaxed_configurations,
            ml_relaxed_slab_result,
            CO_result=ml_gas_references["CO"],
            N2_result=ml_gas_references["N2"],
            H2_result=ml_gas_references["H2"],
            H2O_result=ml_gas_references["H2O"],
        )

    adslab_anomalies_list = [
        job(detect_anomaly)(
            relaxed_result["input_atoms"]["atoms"], relaxed_result["atoms"]
        )
        for relaxed_result in ml_relaxed_configurations
    ]

    top_candidates = filter_sort_select_adslabs(
        adslab_results=ml_relaxed_configurations,
        adslab_anomalies_list=adslab_anomalies_list,
    )

    if num_to_validate_with_DFT == 0:
        return {
            "slab": slab,
            "adslab_ml_relaxed_configurations": top_candidates,
            "adslab_anomalies": adslab_anomalies_list,
        }
    else:
        return {
            "slab": slab,
            "adslab_ml_relaxed_configurations": top_candidates,
            "adslab_anomalies": adslab_anomalies_list,
            "validated_structures": {
                "validated_adslabs": [
                    adslab_validate_job(top_candidates[i]["atoms"])
                    for i in range(num_to_validate_with_DFT)
                ],
                "slab_validated": slab_validate_job(slab.atoms, relax_cell=False),
            },
        }


@job
def reference_adslab_energies(
    adslab_results: list[dict[str, Any]],
    slab_result: RunSchema,
    CO_result: OptSchema | RunSchema,
    H2_result: OptSchema | RunSchema,
    H2O_result: OptSchema | RunSchema,
    N2_result: OptSchema | RunSchema,
):
    adslab_results = copy.deepcopy(adslab_results)
    atomic_energies = {
        "H": H2_result["results"]["energy"] / 2,
        "N": N2_result["results"]["energy"] / 2,
        "O": (H2O_result["results"]["energy"] - H2_result["results"]["energy"]),
        "C": CO_result["results"]["energy"]
        - (H2O_result["results"]["energy"] - H2_result["results"]["energy"]),
    }

    slab_energy = slab_result["results"]["energy"]

    return [
        recursive_dict_merge(
            adslab_result,
            {
                "results": {
                    "atomic_energies": atomic_energies,
                    "slab_energy": slab_energy,
                    "adslab_energy": adslab_result["results"]["energy"],
                    "gas_reactant_energy": sum(
                        [
                            atomic_energies[atom.symbol]
                            for atom in adslab_result["atoms"][
                                adslab_result["atoms"].get_tags() == 2
                            ]  # all adsorbate tagged with tag=2!
                        ]
                    ),
                    "adsorption_energy": adslab_result["results"]["energy"]
                    - slab_energy
                    - sum(
                        [
                            atomic_energies[atom.symbol]
                            for atom in adslab_result["atoms"][
                                adslab_result["atoms"].get_tags() == 2
                            ]  # all adsorbate tagged with tag=2!
                        ]
                    ),
                }
            },
        )
        for adslab_result in adslab_results
    ]


def molecule_pbc(*args, **molecule_kwargs):
    atoms = molecule(*args, **molecule_kwargs)
    atoms.pbc = True
    return atoms


@flow
def standard_ocp_dataset_enumeration(
    bulk_atoms: Atoms,
    adsorbates_kwargs: dict[str, Any],
    multiple_adsorbate_slab_config_kwargs: dict[str, Any],
    ml_relax_job: Job,
    slab_validate_job: Job,
    adslab_validate_job: Job,
    gas_validate_job: Job,
    bulk_relax_job: Job | None,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, dict[str, Any]] | None = None,
    max_miller: int = 1,
    num_to_validate_with_DFT: int = 0,
    reference_ml_energies_to_gas_phase: bool = True,
    relax_bulk: bool = True,
):
    (
        bulk_relax_job_,
        ml_slab_adslab_relax_job_,
        slab_validate_job_,
        adslab_validate_job_,
        gas_validate_job_,
    ) = customize_funcs(
        [
            "bulk_relax_job",
            "ml_slab_adslab_relax_job",
            "slab_validate_job",
            "adslab_validate_job",
            "gas_validate_job",
        ],
        [
            bulk_relax_job,  # type: ignore
            ml_relax_job,
            slab_validate_job,
            adslab_validate_job,
            gas_validate_job,
        ],
        param_swaps=job_params,
        decorators=job_decorators,  # type: ignore
    )

    if relax_bulk:
        bulk_atoms = bulk_relax_job_(bulk_atoms, relax_cell=True)["atoms"]

    slabs = ocp_surface_generator(bulk_atoms=bulk_atoms, max_miller=max_miller)

    if reference_ml_energies_to_gas_phase:
        ml_gas_references = {
            "N2": ml_slab_adslab_relax_job_(
                molecule_pbc("N2", vacuum=10), relax_cell=False
            ),
            "CO": ml_slab_adslab_relax_job_(
                molecule_pbc("CO", vacuum=10), relax_cell=False
            ),
            "H2": ml_slab_adslab_relax_job_(
                molecule_pbc("H2", vacuum=10), relax_cell=False
            ),
            "H2O": ml_slab_adslab_relax_job_(
                molecule_pbc("H2O", vacuum=10), relax_cell=False
            ),
        }
    else:
        ml_gas_references = None

    @flow
    def adsorbML_each_surface(slabs, **kwargs):
        return [
            # We strip the decorator here so it's a bunch of jobs, not subflows. Helpful for prefect!
            strip_decorator(adsorb_ml_pipeline)(slab=slab, **kwargs)
            for slab in slabs
        ]

    return adsorbML_each_surface(
        slabs=slabs,
        adsorbates_kwargs=adsorbates_kwargs,
        multiple_adsorbate_slab_config_kwargs=multiple_adsorbate_slab_config_kwargs,
        ml_slab_adslab_relax_job=ml_slab_adslab_relax_job_,
        slab_validate_job=slab_validate_job_,
        adslab_validate_job=adslab_validate_job_,
        gas_validate_job=gas_validate_job_,
        num_to_validate_with_DFT=num_to_validate_with_DFT,
        ml_gas_references=ml_gas_references,
        reference_ml_energies_to_gas_phase=reference_ml_energies_to_gas_phase,
    )
