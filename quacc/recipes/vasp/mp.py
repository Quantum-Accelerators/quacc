"""
Materials Project-compatible recipes

This set of recipes is meant to be compatible with the Materials Project
Reference: https://doi.org/10.1103/PhysRevMaterials.6.013801
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

    defaults = {"ediffg": -0.05, "xc": "pbesol"}
    flags = defaults | swaps

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
        self.prerelax_kwargs = self.prerelax_kwargs or {}
        self.relax_kwargs = self.relax_kwargs or {}

        # Run the prerelax
        prerelax_results = self.prerelax_electron(atoms, **self.prerelax_kwargs)

        # Update KSPACING arguments
        bandgap = prerelax_results["output"]["bandgap"]
        if bandgap < 1e-4:
            kspacing_swaps = {"kspacing": 0.22, "sigma": 0.2, "ismear": 2, "kpts": None}
        else:
            rmin = 25.22 - 2.87 * bandgap
            kspacing = 2 * np.pi * 1.0265 / (rmin - 1.0183)
            kspacing_swaps = {
                "kspacing": min(kspacing, 0.44),
                "ismear": -5,
                "sigma": 0.05,
                "kpts": None,
            }
        self.relax_kwargs["swaps"] = kspacing_swaps | self.relax_kwargs.get("swaps", {})

        # TODO: Also, copy the WAVECAR from the prerelaxation to the relaxation

        # Run the relax
        return self.relax_electron(prerelax_results["atoms"], **self.relax_kwargs)
