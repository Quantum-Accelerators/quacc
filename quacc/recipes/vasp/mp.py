"""
Materials Project-compatible recipes

This set of recipes is meant to be compatible with the Materials Project
Reference: https://doi.org/10.1038/s41524-022-00881-w
"""
from __future__ import annotations

from dataclasses import dataclass

import covalent as ct
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

        return self.relax_electron(
            self.prerelax_electron(atoms, **prerelax_kwargs)["atoms"], **relax_kwargs
        )
