"""Defect recipes for EMT"""
from __future__ import annotations

from dataclasses import dataclass

import covalent as ct
from ase.atoms import Atoms
from covalent._workflow.electron import Electron
from pymatgen.analysis.defects.generators import AntiSiteGenerator, ChargeInterstitialGenerator, \
    InterstitialGenerator, SubstitutionGenerator, VacancyGenerator, VoronoiInterstitialGenerator

from quacc.recipes.emt.core import relax_job, static_job
from quacc.util.defects import make_defects_from_bulk


@dataclass
class BulkToDefectsFlow:
    """
    Workflow consisting of:

    1. Defect generation

    2. Defect relaxations (optional)

    3. Defect statics (optional)

    Parameters
    ----------
    defect_relax_electron
        Default Electron to use for the relaxation of the defect structures.
    defect_static_electron
        Default Electron to use for the static calculation of the defect structures.
    defect_relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    defect_static_kwargs
        Additional keyword arguments to pass to the static calculation.
    """

    defect_relax_electron: Electron | None = relax_job
    defect_static_electron: Electron | None = static_job
    defect_relax_kwargs: dict | None = None
    defect_static_kwargs: dict | None = None

    def run(
            self,
            atoms: Atoms,
            def_gen: (AntiSiteGenerator | ChargeInterstitialGenerator | InterstitialGenerator | SubstitutionGenerator
                      | VacancyGenerator | VoronoiInterstitialGenerator),
            defectgen_kwargs: dict | None = None,
    ) -> list[dict]:
        """
        Make the workflow.

        Parameters
        ----------
        atoms
            Atoms object for the structure.
        def_gen
            Defect generator
        defectgen_kwargs
            Additional keyword arguments to pass to make_defects_from_bulk()

        Returns
        -------
        list[dict]
            List of dictionary of results from quacc.schemas.ase.summarize_run or quacc.schemas.ase.summarize_opt_run
        """

        self.defect_relax_kwargs = self.defect_relax_kwargs or {}
        self.defect_static_kwargs = self.defect_static_kwargs or {}
        defectgen_kwargs = defectgen_kwargs or {}

        if not self.defect_relax_electron and not self.defect_static_electron:
            raise ValueError(
                "At least one of defect_relax_electron or defect_static_electron must be defined."
            )

        @ct.electron
        @ct.lattice
        def _relax_distributed(defects):
            return [
                self.defect_relax_electron(defect, **self.defect_relax_kwargs)
                for defect in defects
            ]

        @ct.electron
        @ct.lattice
        def _static_distributed(defects):
            return [
                self.defect_static_electron(defect, **self.defect_static_kwargs)
                for defect in defects
            ]

        @ct.electron
        @ct.lattice
        def _relax_and_static_distributed(defects):
            return [
                self.defect_static_electron(
                    self.defect_relax_electron(defect, **self.defect_relax_kwargs)["atoms"],
                    **self.defect_static_kwargs,
                )
                for defect in defects
            ]

        defects = ct.electron(make_defects_from_bulk)(atoms, def_gen, **defectgen_kwargs)

        if self.defect_relax_electron and self.defect_static_electron:
            return _relax_and_static_distributed(defects)
        elif self.defect_relax_electron:
            return _relax_distributed(defects)
        elif self.defect_static_electron:
            return _static_distributed(defects)
