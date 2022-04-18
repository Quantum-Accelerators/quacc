"""Core recipes for DFTB+"""
from __future__ import annotations
from dataclasses import dataclass, field
from shutil import which
from typing import Any, Dict, List

from ase.atoms import Atoms
from ase.calculators.dftb import Dftb
from jobflow import Maker, job
from monty.dev import requires

from quacc.schemas.calc import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc

DFTBPLUS_EXISTS = bool(which("dftb+"))
GEOM_FILE = "geo_end.gen"


@dataclass
class StaticJob(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    method
        Method to use. Accepts 'DFTB', 'GFN1-xTB', and 'GFN2-xTB'.
    kpts
        k-point grid to use. Defaults to None for molecules and
        (1, 1, 1) for solids.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "DFTB-Static"
    method: str = "GFN2-xTB"
    kpts: tuple | List[tuple] | Dict[str, Any] = None
    swaps: Dict[str, Any] = field(default_factory=dict)

    @job
    @requires(
        DFTBPLUS_EXISTS,
        "DFTB+ must be installed. Try conda install -c conda-forge dftbplus",
    )
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object

        Returns
        -------
        Dict
            Summary of the run.
        """
        defaults = {}
        if "xtb" in self.method.lower():
            defaults["Hamiltonian_"] = "xTB"
            if "gfn2-xtb" in self.method.lower():
                defaults["Hamiltonian_Method"] = "GFN2-xTB"
            elif "gfn1-xtb" in self.method.lower():
                defaults["Hamiltonian_Method"] = "GFN1-xTB"

        if self.kpts:
            defaults["kpts"] = self.kpts
        elif atoms.pbc.any():
            defaults["kpts"] = (1, 1, 1)
        else:
            defaults["kpts"] = None

        flags = merge_dicts(
            defaults, self.swaps, remove_none=True, auto_lowercase=False
        )

        atoms.calc = Dftb(**flags)
        new_atoms = run_calc(atoms, geom_file=GEOM_FILE)
        summary = summarize_run(
            new_atoms, input_atoms=atoms, additional_fields={"name": self.name}
        )

        return summary


@dataclass
class RelaxJob(Maker):
    """
    Class to carry out a structure relaxation.

    Parameters
    ----------
    name
        Name of the job.
    method
        Method to use. Accepts 'DFTB', 'GFN1-xTB', and 'GFN2-xTB'.
    kpts
        k-point grid to use. Defaults to None for molecules and
        (1, 1, 1) for solids.
    lattice_opt
        Whether to relax the unit cell shape/volume in addition to
        the positions.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "DFTB-Relax"
    method: str = "GFN2-xTB"
    kpts: tuple | List[tuple] | Dict[str, Any] = None
    lattice_opt: bool = False
    swaps: Dict[str, Any] = field(default_factory=dict)

    @job
    @requires(
        DFTBPLUS_EXISTS,
        "DFTB+ must be installed. Try conda install -c conda-forge dftbplus",
    )
    def make(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Make the run.

        Parameters
        ----------
        atoms
            .Atoms object

        Returns
        -------
        Dict
            Summary of the run.
        """
        defaults = {
            "Driver_": "GeometryOptimization",
            "Driver_LatticeOpt": "Yes"
            if self.lattice_opt and True in atoms.pbc
            else "No",
            "Driver_MaxSteps": 1000,
        }

        if "xtb" in self.method.lower():
            defaults["Hamiltonian_"] = "xTB"
            if "gfn2-xtb" in self.method.lower():
                defaults["Hamiltonian_Method"] = "GFN2-xTB"
            elif "gfn1-xtb" in self.method.lower():
                defaults["Hamiltonian_Method"] = "GFN1-xTB"

        if self.kpts:
            defaults["kpts"] = self.kpts
        elif atoms.pbc.any():
            defaults["kpts"] = (1, 1, 1)
        else:
            defaults["kpts"] = None

        flags = merge_dicts(
            defaults, self.swaps, remove_none=True, auto_lowercase=False
        )

        atoms.calc = Dftb(**flags)
        new_atoms = run_calc(atoms, geom_file=GEOM_FILE)
        summary = summarize_run(
            new_atoms, input_atoms=atoms, additional_fields={"name": self.name}
        )

        return summary
