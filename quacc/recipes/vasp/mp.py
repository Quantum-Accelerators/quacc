"""
Materials Project-compatible recipes

This set of recipes is meant to be compatible with the Materials Project
Reference: https://doi.org/10.1038/s41524-022-00881-w
"""
from __future__ import annotations

from dataclasses import dataclass

import covalent as ct
import numpy as np
from ase import Atoms
from covalent._workflow.electron import Electron

from quacc.calculators.vasp import Vasp
from quacc.schemas.vasp import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import merge_dicts


@ct.electron
def mp_prerelax_job(
    atoms: Atoms, preset: str | None = "MPScanSet", swaps: dict = None
) -> dict:
    """
    Function to pre-relax a structure with Materials Project settings.
    By default, this uses a PBEsol pre-relax step.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.vasp.summarize_run
    """
    swaps = swaps or {}

    defaults = {"xc": "pbesol", "ediffg": -0.05}
    flags = merge_dicts(defaults, swaps)

    calc = Vasp(atoms, preset=preset, **flags)
    atoms.calc = calc
    atoms = run_calc(atoms)

    return summarize_run(atoms, additional_fields={"name": "MP-Prerelax"})


@ct.electron
def mp_relax_job(
    atoms: Atoms, preset: str | None = "MPScanSet", swaps: dict = None
) -> dict:
    """
    Function to relax a structure with Materials Project settings.
    By default, this uses an r2SCAN relax step.

    Parameters
    ----------
    atoms
        Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from quacc.schemas.vasp.summarize_run
    """
    swaps = swaps or {}

    calc = Vasp(atoms, preset=preset, **swaps)
    atoms.calc = calc
    atoms = run_calc(atoms)

    return summarize_run(atoms, additional_fields={"name": "MP-Relax"})


@dataclass
class MPRelaxFlow:
    """
    Workflow consisting of:

    1. MP-compatible pre-relax

    2. MP-compatible relax

    Parameters
    ----------
    prerelax_electron
        Default to use for the pre-relaxation.
    relax_electron
        Default to use for the relaxation.
    prerelax_kwargs
        Additional keyword arguments to pass to the pre-relaxation calculation.
    relax_kwargs
        Additional keyword arguments to pass to the relaxation calculation.
    """

    prerelax_electron: Electron | None = mp_prerelax_job
    relax_electron: Electron | None = mp_relax_job
    prerelax_kwargs: dict = None
    relax_kwargs: dict = None

    def run(self, atoms: Atoms) -> dict:
        """
        Run the workflow.

        Parameters
        ----------
        atoms
            Atoms object for the structure.

        Returns
        -------
        dict
            Dictionary results from quacc.schemas.vasp.summarize_run
        """
        prerelax_kwargs = self.prerelax_kwargs or {}
        relax_kwargs = self.relax_kwargs or {}

        # TODO: Also, copy the WAVECAR
        prerelax_results = self.prerelax_electron(atoms, **prerelax_kwargs)["atoms"]
        self._set_kspacing_swaps(prerelax_results["output"]["bandgap"])

        return self.relax_electron(prerelax_results["atoms"], **relax_kwargs)

    def _set_kspacing_swaps(self, bandgap: float) -> dict:
        """
        Function to calculate KSPACING and related parameters for a given bandgap.

        Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801

        Parameters
        ----------
        bandgap
            Bandgap of the structure in eV.

        Returns
        -------
        dict
            Dictionary of swaps.
        """

        if bandgap < 1e-4:
            kspacing_swaps = {"kspacing": 0.22, "sigma": 0.2, "ismear": 2, "kpts": None}
        else:
            rmin = 25.22 - 2.87 * bandgap
            kspacing = 2 * np.pi * 1.0265 / (rmin - 1.0183)
            kspacing_swaps = {
                "kspacing": kspacing if 0.22 < kspacing < 0.44 else 0.44,
                "ismear": -5,
                "sigma": 0.05,
                "kpts": None,
            }

        self.relax_kwargs["swaps"] = merge_dicts(
            kspacing_swaps, self.relax_kwargs.get("swaps", {})
        )
