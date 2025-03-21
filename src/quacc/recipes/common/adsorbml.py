from __future__ import annotations

import copy
import logging
from typing import TYPE_CHECKING, Any, Literal, TypedDict

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

class AdsorbatesKwargs(TypedDict):
    adsorbate_type: str
    position: list[float]
    orientation: list[float]

class MultipleAdsorbateSlabConfigKwargs(TypedDict):
    num_adsorbates: int
    adsorbate_distance: float

class MoleculeResults(TypedDict):
    N2: RunSchema
    CO: RunSchema
    H2: RunSchema
    H2O: RunSchema

class AtomicReferenceEnergies(TypedDict):
    H: float
    N: float
    O: float
    C: float

@job
def ocp_surface_generator(bulk_atoms: Atoms, max_miller: int = 1) -> list[Slab]:
    """
    Generate surface slabs from bulk atoms.

    Parameters
    ----------
    bulk_atoms : Atoms
        The bulk atomic structure.
    max_miller : int, optional
        Maximum Miller index, by default 1.

    Returns
    -------
    list[Slab]
        List of generated surface slabs.
    """
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
        """
        Initialize a CustomSlab object.

        Parameters
        ----------
        slab_atoms : Atoms
            The slab atomic structure.
        bulk : Atoms, optional
            The bulk atomic structure, by default None.
        millers : tuple, optional
            Miller indices, by default None.
        shift : float, optional
            Shift value, by default None.
        top : bool, optional
            Top value, by default None.
        min_ab : float, optional
            Minimum a and b lattice parameters, by default 0.8.
        """
        self.bulk = bulk
        self.slab_atoms = slab_atoms
        self.millers = millers
        self.shift = shift
        self.top = top

        assert np.linalg.norm(self.slab_atoms.cell[0]) >= min_ab, "Slab not tiled"
        assert np.linalg.norm(self.slab_atoms.cell[1]) >= min_ab, "Slab not tiled"
        assert self.has_surface_tagged(), "Slab not tagged"
        assert len(self.slab_atoms.constraints) > 0, "Sub-surface atoms not constrained"


@job
def ocp_adslab_generator(
    slab: Slab | Atoms,
    adsorbates_kwargs: list[AdsorbatesKwargs] | None = None,
    multiple_adsorbate_slab_config_kwargs: MultipleAdsorbateSlabConfigKwargs | None = None,
) -> list[Atoms]:
    """
    Generate adsorbate-slab configurations.

    Parameters
    ----------
    slab : Slab | Atoms
        The slab structure.
    adsorbates_kwargs : list[AdsorbatesKwargs], optional
        List of keyword arguments for generating adsorbates, by default None.
    multiple_adsorbate_slab_config_kwargs : MultipleAdsorbateSlabConfigKwargs, optional
        Keyword arguments for generating multiple adsorbate-slab configurations, by default None.

    Returns
    -------
    list[Atoms]
        List of generated adsorbate-slab configurations.
    """
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
    adsorbates_kwargs: AdsorbatesKwargs,
    multiple_adsorbate_slab_config_kwargs: MultipleAdsorbateSlabConfigKwargs | None = None,
) -> list[dict[str, Any]]:
    """
    Find adsorbate-slab configurations for each slab.

    Parameters
    ----------
    slabs : list[Slab]
        List of slabs.
    adsorbates_kwargs : AdsorbatesKwargs
        Keyword arguments for generating adsorbates.
    multiple_adsorbate_slab_config_kwargs : MultipleAdsorbateSlabConfigKwargs, optional
        Keyword arguments for generating multiple adsorbate-slab configurations, by default None.

    Returns
    -------
    list[dict[str, Any]]
        List of dictionaries containing slabs and their corresponding adsorbate-slab configurations.
    """
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
    """
    Detect anomalies between initial and final atomic structures.

    Parameters
    ----------
    initial_atoms : Atoms
        Initial atomic structure.
    final_atoms : Atoms
        Final atomic structure.

    Returns
    -------
    list[Literal["adsorbate_dissociated", "adsorbate_desorbed", "surface_changed", "adsorbate_intercalated"]]
        List of detected anomalies.
    """
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
) -> list[RunSchema]:
    """
    Filter, sort, and select adsorbate-slab configurations based on anomalies and energy.

    Parameters
    ----------
    adslab_results : list[RunSchema]
        List of adsorbate-slab results.
    adslab_anomalies_list : list[list[str]]
        List of detected anomalies for each adsorbate-slab configuration.

    Returns
    -------
    list[RunSchema]
        Sorted list of adsorbate-slab configurations without anomalies.
    """
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
    adsorbates_kwargs: AdsorbatesKwargs,
    multiple_adsorbate_slab_config_kwargs: MultipleAdsorbateSlabConfigKwargs,
    ml_slab_adslab_relax_job: Job,
    slab_validate_job: Job,
    adslab_validate_job: Job,
    gas_validate_job: Job,
    num_to_validate_with_DFT: int = 0,
    reference_ml_energies_to_gas_phase: bool = False,
    molecule_results: MoleculeResults | None = None,
    atomic_reference_energies: AtomicReferenceEnergies | None = None,
) -> dict[str, Any]:
    """
    Run a machine learning-based pipeline for adsorbate-slab systems.

    Parameters
    ----------
    slab : Slab
        The slab structure to which adsorbates will be added.
    adsorbates_kwargs : AdsorbatesKwargs
        Keyword arguments for generating adsorbate configurations.
    multiple_adsorbate_slab_config_kwargs : MultipleAdsorbateSlabConfigKwargs
        Keyword arguments for generating multiple adsorbate-slab configurations.
    ml_slab_adslab_relax_job : Job
        Job for relaxing slab and adsorbate-slab configurations using ML.
    slab_validate_job : Job
        Job for validating the slab structure.
    adslab_validate_job : Job
        Job for validating the adsorbate-slab structures.
    gas_validate_job : Job
        Job for validating gas phase structures.
    num_to_validate_with_DFT : int, optional
        Number of top configurations to validate with DFT, by default 0.
    reference_ml_energies_to_gas_phase : bool, optional
        Whether to reference ML energies to gas phase, by default False.
    molecule_results : MoleculeResults, optional
        Precomputed molecule results for referencing, by default None.
    atomic_reference_energies : AtomicReferenceEnergies, optional
        Atomic reference energies for referencing, by default None.

    Returns
    -------
    dict
        Dictionary containing the slab, ML-relaxed adsorbate-slab configurations,
        detected anomalies, and optionally DFT-validated structures.
    """
    unrelaxed_adslab_configurations = ocp_adslab_generator(
        slab, adsorbates_kwargs, multiple_adsorbate_slab_config_kwargs
    )

    ml_relaxed_configurations = [
        ml_slab_adslab_relax_job(adslab_configuration)
        for adslab_configuration in unrelaxed_adslab_configurations
    ]

    ml_relaxed_slab_result = ml_slab_adslab_relax_job(slab.atoms)

    if reference_ml_energies_to_gas_phase:
        if atomic_reference_energies is None and molecule_results is None:
            molecule_results = generate_molecule_reference_results(
                ml_slab_adslab_relax_job
            )

        ml_relaxed_configurations = reference_adslab_energies(
            ml_relaxed_configurations,
            ml_relaxed_slab_result,
            atomic_energies=atomic_reference_energies,
            molecule_results=molecule_results,
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
    atomic_energies: AtomicReferenceEnergies | None,
    molecule_results: MoleculeResults | None,
) -> list[dict[str, Any]]:
    """
    Reference adsorbate-slab energies to atomic and slab energies.

    Parameters
    ----------
    adslab_results : list[dict[str, Any]]
        List of adsorbate-slab results.
    slab_result : RunSchema
        Result of the slab calculation.
    atomic_energies : AtomicReferenceEnergies | None
        Dictionary of atomic energies.
    molecule_results : MoleculeResults | None
        Dictionary of molecule results.

    Returns
    -------
    list[dict[str, Any]]
        List of adsorbate-slab results with referenced energies.
    """
    adslab_results = copy.deepcopy(adslab_results)
    if atomic_energies is None and molecule_results is not None:
        atomic_energies = {
            "H": molecule_results["H2"]["results"]["energy"] / 2,
            "N": molecule_results["N2"]["results"]["energy"] / 2,
            "O": (
                molecule_results["H2O"]["results"]["energy"]
                - molecule_results["H2"]["results"]["energy"]
            ),
            "C": molecule_results["CO"]["results"]["energy"]
            - (
                molecule_results["H2O"]["results"]["energy"]
                - molecule_results["H2"]["results"]["energy"]
            ),
        }
    else:
        raise Exception(
            "Missing atomic energies and gas phase energies; unable to continue!"
        )

    slab_energy = slab_result["results"]["energy"]

    return [
        recursive_dict_merge(
            adslab_result,
            {
                "results": {
                    "referenced_adsorption_energy": {
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
                }
            },
        )
        for adslab_result in adslab_results
    ]


def molecule_pbc(*args: Any, **molecule_kwargs: Any) -> Atoms:
    """
    Create a molecule with periodic boundary conditions.

    Parameters
    ----------
    *args : Any
        Positional arguments for the molecule function.
    **molecule_kwargs : Any
        Keyword arguments for the molecule function.

    Returns
    -------
    Atoms
        Atomic structure with periodic boundary conditions.
    """
    atoms = molecule(*args, **molecule_kwargs)
    atoms.pbc = True
    return atoms


def generate_molecule_reference_results(relax_job: Job) -> MoleculeResults:
    """
    Generate reference results for molecules.

    Parameters
    ----------
    relax_job : Job
        Job for relaxing molecular structures.

    Returns
    -------
    MoleculeResults
        Dictionary of reference results for molecules.
    """
    return {
        "N2": relax_job(molecule_pbc("N2", vacuum=10), relax_cell=False),
        "CO": relax_job(molecule_pbc("CO", vacuum=10), relax_cell=False),
        "H2": relax_job(molecule_pbc("H2", vacuum=10), relax_cell=False),
        "H2O": relax_job(molecule_pbc("H2O", vacuum=10), relax_cell=False),
    }


@flow
def bulk_to_surfaces_to_adsorbml(
    bulk_atoms: Atoms,
    adsorbates_kwargs: AdsorbatesKwargs,
    multiple_adsorbate_slab_config_kwargs: MultipleAdsorbateSlabConfigKwargs,
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
) -> list[dict[str, Any]]:
    """
    Run a pipeline from bulk atoms to adsorbate-slab configurations using machine learning.

    Parameters
    ----------
    bulk_atoms : Atoms
        The bulk atomic structure.
    adsorbates_kwargs : AdsorbatesKwargs
        Keyword arguments for generating adsorbate configurations.
    multiple_adsorbate_slab_config_kwargs : MultipleAdsorbateSlabConfigKwargs
        Keyword arguments for generating multiple adsorbate-slab configurations.
    ml_relax_job : Job
        Job for relaxing slab and adsorbate-slab configurations using ML.
    slab_validate_job : Job
        Job for validating the slab structure.
    adslab_validate_job : Job
        Job for validating the adsorbate-slab structures.
    gas_validate_job : Job
        Job for validating gas phase structures.
    bulk_relax_job : Job | None
        Job for relaxing the bulk structure, by default None.
    job_params : dict[str, dict[str, Any]], optional
        Parameters for customizing jobs, by default None.
    job_decorators : dict[str, dict[str, Any]], optional
        Decorators for customizing jobs, by default None.
    max_miller : int, optional
        Maximum Miller index, by default 1.
    num_to_validate_with_DFT : int, optional
        Number of top configurations to validate with DFT, by default 0.
    reference_ml_energies_to_gas_phase : bool, optional
        Whether to reference ML energies to gas phase, by default True.
    relax_bulk : bool, optional
        Whether to relax the bulk structure, by default True.

    Returns
    -------
    list[dict[str, Any]]
        List of dictionaries containing the results of the pipeline for each slab.
    """
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
        molecule_results = generate_molecule_reference_results(
            ml_slab_adslab_relax_job_
        )
    else:
        molecule_results = None

    @flow
    def adsorbML_each_surface(slabs: list[Slab], **kwargs: Any) -> list[dict[str, Any]]:
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
        molecule_results=molecule_results,
        reference_ml_energies_to_gas_phase=reference_ml_energies_to_gas_phase,
    )
