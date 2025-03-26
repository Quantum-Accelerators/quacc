from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

from quacc import Job, flow
from quacc.recipes.common.adsorbml import (
    adsorb_ml_pipeline as adsorb_ml_pipeline_,
    bulk_to_surfaces_to_adsorbml as bulk_to_surfaces_to_adsorbml_,
)
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from fairchem.data.oc.core.slab import Slab

    from quacc.types import (
        AdsorbMLSchema,
        AtomicReferenceEnergies,
        MoleculeReferenceResults,
    )



logger = logging.getLogger(__name__)


@flow
def adsorb_ml_pipeline(
    slab: Slab,
    adsorbates_kwargs: dict[str, Any],
    multiple_adsorbate_slab_config_kwargs: dict[str, Any],
    ml_slab_adslab_relax_job: Job,
    slab_validate_job: Job,
    adslab_validate_job: Job,
    gas_validate_job: Job,
    num_to_validate_with_DFT: int = 0,
    reference_ml_energies_to_gas_phase: bool = False,
    molecule_results: MoleculeReferenceResults | None = None,
    atomic_reference_energies: AtomicReferenceEnergies | None = None,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, dict[str, Any]] | None = None,
) -> AdsorbMLSchema:
    """
    Run a machine learning-based pipeline for adsorbate-slab systems.

    1. Relax slab using ML
    2. Generate trial adsorbate-slab configurations for the relaxed slab
    3. Relax adsorbate-slab configurations using ML
    4. Validate slab and adsorbate-slab configurations (check for anomalies like dissociations))
    5. Reference the energies to gas phase if needed (eg using a total energy ML model)
    6. Optionally validate top K configurations with DFT single-points or relaxations

    Parameters
    ----------
    slab : Slab
        The slab structure to which adsorbates will be added.
    adsorbates_kwargs : dict[str, Any]
        Keyword arguments for generating adsorbate configurations.
    multiple_adsorbate_slab_config_kwargs : dict[str, Any]
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
    molecule_results : MoleculeReferenceResults, optional
        Precomputed molecule results for referencing, by default None.
    atomic_reference_energies : AtomicReferenceEnergies, optional
        Atomic reference energies for referencing, by default None.

    Returns
    -------
    dict
        Dictionary containing the slab, ML-relaxed adsorbate-slab configurations,
        detected anomalies, and optionally DFT-validated structures.
    """

    (
        ml_slab_adslab_relax_job_,
        slab_validate_job_,
        adslab_validate_job_,
        gas_validate_job_,
    ) = customize_funcs(
        [
            "ml_slab_adslab_relax_job",
            "slab_validate_job",
            "adslab_validate_job",
            "gas_validate_job",
        ],
        [
            ml_slab_adslab_relax_job,
            slab_validate_job,
            adslab_validate_job,
            gas_validate_job,
        ],
        param_swaps=job_params,
        decorators=job_decorators,  # type: ignore
    )

    return adsorb_ml_pipeline_(
        slab,
        adsorbates_kwargs,
        multiple_adsorbate_slab_config_kwargs,
        ml_slab_adslab_relax_job_,
        slab_validate_job_,
        adslab_validate_job_,
        gas_validate_job_,
        num_to_validate_with_DFT,
        reference_ml_energies_to_gas_phase,
        molecule_results,
        atomic_reference_energies,
    )


@flow
def bulk_to_surfaces_to_adsorbml(
    bulk_atoms: Atoms,
    adsorbates_kwargs: dict[str,Any],
    multiple_adsorbate_slab_config_kwargs: dict[str, Any],
    ml_relax_job: Job,
    slab_validate_job: Job,
    adslab_validate_job: Job,
    gas_validate_job: Job,
    max_miller: int = 1,
    bulk_relax_job: Job | None = None,
    job_params: dict[str, dict[str, Any]] | None = None,
    job_decorators: dict[str, dict[str, Any]] | None = None,
    num_to_validate_with_DFT: int = 0,
    reference_ml_energies_to_gas_phase: bool = True,
    relax_bulk: bool = True,
    atomic_reference_energies: dict[str, float] | None = None,
) -> list[dict[str, Any]]:
    """
    Run a pipeline from bulk atoms to adsorbate-slab configurations using machine learning!
    For full details, see the AdsorbML paper (https://arxiv.org/abs/2211.16486,
                                     https://www.nature.com/articles/s41524-023-01121-5).

    1. Relax bulk structure if desired
    2. Generate surface slabs from bulk atoms
    3. Generate gas phase reference energies if needed

    For each slab generated in (3):
        1. Relax slab using ML
        2. Generate trial adsorbate-slab configurations for the relaxed slab
        3. Relax adsorbate-slab configurations using ML
        4. Validate slab and adsorbate-slab configurations (check for anomalies like dissociations))
        5. Reference the energies to gas phase if needed (eg using a total energy ML model)
        6. Optionally validate top K configurations with DFT single-points or relaxations

    Parameters
    ----------
    bulk_atoms : Atoms
        The bulk atomic structure.
    adsorbates_kwargs : AdsorbatesKwargs
        Keyword arguments for generating adsorbate configurations.
    multiple_adsorbate_slab_config_kwargs : dict[str, Any]
        Keyword arguments for generating multiple adsorbate-slab configurations.
    ml_relax_job : Job
        Job for relaxing slab and adsorbate-slab configurations using ML.
    slab_validate_job : Job
        Job for validating the slab structure.
    adslab_validate_job : Job
        Job for validating the adsorbate-slab structures.
    gas_validate_job : Job
        Job for validating gas phase structures.
    max_miller : int, optional
        Maximum Miller index, by default 1.
    bulk_relax_job : Job | None, optional
        Job for relaxing the bulk structure, by default None.
    job_params : dict[str, dict[str, Any]] | None, optional
        Parameters for customizing jobs, by default None.
    job_decorators : dict[str, dict[str, Any]] | None, optional
        Decorators for customizing jobs, by default None.
    num_to_validate_with_DFT : int, optional
        Number of top configurations to validate with DFT, by default 0.
    reference_ml_energies_to_gas_phase : bool, optional
        Whether to reference ML energies to gas phase, by default True.
    relax_bulk : bool, optional
        Whether to relax the bulk structure, by default True.
    atomic_reference_energies : dict[str, float] | None, optional
        Atomic reference energies for referencing if known ahead of time, by default None.

    Returns
    -------
    list[AdsorbMLSchema]
        List of AdsorbML results for each slab
    """
    (
        bulk_relax_job_,
        ml_relax_job_,
        slab_validate_job_,
        adslab_validate_job_,
        gas_validate_job_,
    ) = customize_funcs(
        [
            "bulk_relax_job",
            "ml_relax_job",
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
    return bulk_to_surfaces_to_adsorbml_(
        bulk_atoms,
        adsorbates_kwargs,
        multiple_adsorbate_slab_config_kwargs,
        ml_relax_job_,
        slab_validate_job_,
        adslab_validate_job_,
        gas_validate_job_,
        max_miller,
        bulk_relax_job_,
        num_to_validate_with_DFT,
        reference_ml_energies_to_gas_phase,
        relax_bulk,
        atomic_reference_energies,
    )
