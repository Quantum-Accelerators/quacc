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
    ) = VacancyGenerator,  # NOTE: I added a default. VacancyGenerator seemed like it might be the most popular.
    charge_state: int = 0,  # NOTE: I changed this from int | None = None because None does not work
    make_defects_kwargs: dict | None = None,
    defect_relax: ct.electron | None = relax_job,
    defect_static: ct.electron | None = static_job,
    defect_relax_kwargs: dict | None = None,
    defect_static_kwargs: dict | None = None,
) -> list[RunSchema | OptSchema]:
    """
    Workflow consisting of:

    1. Defect generation

    2. Defect relaxations

    3. Defect statics (optional)

    Parameters
    ----------
    atoms
        Atoms object for the structure.
    defectgen
        Defect generator
    charge_state
        Charge state of the defect
    make_defects_kwargs
        Keyword arguments to pass to the make_defects_from_bulk
    defect_relax
        Default Electron to use for the relaxation of the defect structures.
    defect_static
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
    defect_relax_kwargs = defect_relax_kwargs or {"relax_cell": False}
    defect_static_kwargs = defect_static_kwargs or {}
    make_defects_kwargs = make_defects_kwargs or {}

    @ct.electron
    def _make_defects(atoms):
        atoms = atoms if isinstance(atoms, Atoms) else atoms["atoms"]
        return make_defects_from_bulk(
            atoms, defectgen=defectgen, charge_state=charge_state, **make_defects_kwargs
        )

    @ct.electron
    @ct.lattice
    def _relax_distributed(defects):
        return [defect_relax(defect, **defect_relax_kwargs) for defect in defects]

    @ct.electron
    @ct.lattice
    def _relax_and_static_distributed(defects):
        return [
            defect_static(
                defect_relax(defect, **defect_relax_kwargs),
                **defect_static_kwargs,
            )
            for defect in defects
        ]

    defects = _make_defects(atoms)

    if defect_static is None:
        return _relax_distributed(defects)

    return _relax_and_static_distributed(defects)
