"""Defect recipes for EMT"""
from __future__ import annotations

import covalent as ct
from ase.atoms import Atoms
from pymatgen.analysis.defects.generators import (
    AntiSiteGenerator,
    ChargeInterstitialGenerator,
    InterstitialGenerator,
    SubstitutionGenerator,
    VacancyGenerator,
    VoronoiInterstitialGenerator,
)

from quacc.recipes.emt.core import relax_job, static_job
from quacc.schemas.ase import OptSchema, RunSchema
from quacc.util.defects import make_defects_from_bulk


def bulk_to_defects_flow(
    atoms: Atoms | dict,
    defectgen: (
        AntiSiteGenerator
        | ChargeInterstitialGenerator
        | InterstitialGenerator
        | SubstitutionGenerator
        | VacancyGenerator
        | VoronoiInterstitialGenerator
    ),
    charge_state: int | None = None,
    defectgen_kwargs: dict | None = None,
    defect_relax_electron: ct.electron | None = relax_job,
    defect_static_electron: ct.electron | None = static_job,
    defect_relax_kwargs: dict | None = None,
    defect_static_kwargs: dict | None = None,
) -> list[RunSchema | OptSchema]:
    """
    Workflow consisting of:

    1. Defect generation

    2. Defect relaxations (optional)

    3. Defect statics (optional)

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    defectgen
        Defect generator
    charge_state
        Charge state of the defect
    defectgen_kwargs
        Keyword arguments to pass to the pymatgen.analysis.defects.generators.get_defects() method
    defect_relax_electron
        Default Electron to use for the relaxation of the defect structures.
    defect_static_electron
        Default Electron to use for the static calculation of the defect structures.
    defect_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    defect_static_kwargs
        Additional keyword arguments to pass to the static calculation.

    Returns
    -------
    list[dict]
        List of dictionary of results from quacc.schemas.ase.summarize_run or quacc.schemas.ase.summarize_opt_run
    """
    atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
    defect_relax_kwargs = defect_relax_kwargs or {}
    defect_static_kwargs = defect_static_kwargs or {}
    defectgen_kwargs = defectgen_kwargs or {}

    if not defect_relax_electron and not defect_static_electron:
        raise ValueError(
            "At least one of defect_relax_electron or defect_static_electron must be defined."
        )

    @ct.electron
    @ct.lattice
    def _relax_distributed(defects):
        return [
            defect_relax_electron(defect, **defect_relax_kwargs) for defect in defects
        ]

    @ct.electron
    @ct.lattice
    def _static_distributed(defects):
        return [
            defect_static_electron(defect, **defect_static_kwargs) for defect in defects
        ]

    @ct.electron
    @ct.lattice
    def _relax_and_static_distributed(defects):
        return [
            defect_static_electron(
                defect_relax_electron(defect, **defect_relax_kwargs)["atoms"],
                **defect_static_kwargs,
            )
            for defect in defects
        ]

    defects = ct.electron(make_defects_from_bulk)(
        atoms, defectgen, charge_state, **defectgen_kwargs
    )

    if defect_relax_electron and defect_static_electron:
        return _relax_and_static_distributed(defects)
    elif defect_relax_electron:
        return _relax_distributed(defects)
    elif defect_static_electron:
        return _static_distributed(defects)
