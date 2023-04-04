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
from quacc.util.calc import _check_logfile, run_calc

DFTBPLUS_EXISTS = bool(which("dftb+"))
LOG_FILE = "dftb.out"
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
        defaults = {
            "Hamiltonian_": "xTB" if "xtb" in self.method.lower() else "DFTB",
            "Hamiltonian_Method": self.method if "xtb" in self.method.lower() else None,
            "kpts": self.kpts if self.kpts else (1, 1, 1) if atoms.pbc.any() else None,
        }
        flags = merge_dicts(
            defaults, self.swaps, remove_none=True, auto_lowercase=False
        )

        atoms.calc = Dftb(**flags)
        new_atoms = run_calc(atoms, geom_file=GEOM_FILE)
        scc_check = _check_logfile(LOG_FILE, "SCC is NOT converged")
        if scc_check:
            raise ValueError("SCC is not converged")
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
            "Hamiltonian_": "xTB" if "xtb" in self.method.lower() else "DFTB",
            "Hamiltonian_Method": self.method if "xtb" in self.method.lower() else None,
            "kpts": self.kpts if self.kpts else (1, 1, 1) if atoms.pbc.any() else None,
            "Driver_": "GeometryOptimization",
            "Driver_LatticeOpt": "Yes" if self.lattice_opt else "No",
            "Driver_AppendGeometries": "Yes",
            "Driver_MaxSteps": 2000,
        }
        flags = merge_dicts(
            defaults, self.swaps, remove_none=True, auto_lowercase=False
        )

        atoms.calc = Dftb(**flags)
        new_atoms = run_calc(atoms, geom_file=GEOM_FILE)
        geom_check = _check_logfile(LOG_FILE, "Geometry converged")
        if not geom_check:
            raise ValueError("Geometry did not converge")
        summary = summarize_run(
            new_atoms, input_atoms=atoms, additional_fields={"name": self.name}
        )

        return summary
