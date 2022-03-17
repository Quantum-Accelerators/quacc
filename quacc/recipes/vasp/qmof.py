from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from ase.optimize import BFGSLineSearch
from jobflow import Maker, job

from quacc.calculators.vasp import SmartVasp
from quacc.schemas.calc import summarize_run as summarize_ase_run
from quacc.schemas.vasp import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc

# This set of recipes is meant to be compatible with the QMOF Database workflow.
# Reference: https://doi.org/10.1016/j.matt.2021.02.015


@dataclass
class QMOFJob(Maker):
    """
    Class to relax a structure in a multi-step process for increased
    computational efficiency. This is all done in a single compute job.
    Settings are such that they are compatible with the QMOF Database.

    1. A "pre-relaxation" with BFGSLineSearch to resolve very high forces.
    2. Position relaxation with default ENCUT and coarse k-point grid.
    3. Optional: volume relaxation with coarse k-point grid.
    4. Double relaxation using production-quality settings.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use. Applies for all jobs.
    volume_relax
        True if a volume relaxation should be performed.
        False if only the positions should be updated.
    swaps
        Dictionary of custom kwargs for the calculator. Applies for all jobs.
    """

    name: str = "QMOF"
    preset: str = "QMOFSet"
    volume_relax: bool = True
    swaps: Dict[str, Any] = None

    @job
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
        swaps = self.swaps or {}

        # 1. Pre-relaxation
        summary1 = prerelax(atoms, self.preset, swaps, fmax=5.0)
        atoms = summary1["atoms"]

        # 2. Position relaxation (loose)
        summary2 = loose_relax_positions(atoms, self.preset, swaps)
        atoms = summary2["atoms"]

        # 3. Optional: Volume relaxation (loose)
        if self.volume_relax:
            summary3 = loose_relax_volume(atoms, self.preset, swaps)
            atoms = summary3["atoms"]

        # 4. Double Relaxation
        # This is done for two reasons: a) because it can resolve repadding
        # issues when dV is large; b) because we can use LREAL = Auto for the
        # first relaxation and the default LREAL for the second.
        summary4 = double_relax(
            atoms, self.preset, swaps, volume_relax=self.volume_relax
        )
        atoms = summary4["relax2"]["atoms"]

        # 5. Static Calculation
        summary5 = static(atoms, self.preset, swaps)

        return {
            "prerelax-lowacc": summary1,
            "position-relax-lowacc": summary2,
            "volume-relax-lowacc": summary3 if self.volume_relax else None,
            "double-relax": summary4,
            "static": summary5,
        }


def prerelax(
    atoms: Atoms,
    preset: str = "QMOFSet",
    swaps: Dict[str, Any] = None,
    fmax: float = 5.0,
) -> Dict[str, Any]:
    """
    A "pre-relaxation" with BFGSLineSearch to resolve very high forces.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    fmax
        Maximum force in eV/A.

    Returns
    -------
    Dict
        Summary of the run.
    """
    swaps = swaps or {}
    defaults = {
        "auto_kpts": {"grid_density": 100},
        "ediff": 1e-4,
        "encut": None,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nelm": 225,
        "nsw": 0,
    }
    flags = merge_dicts(defaults, swaps)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    dyn = BFGSLineSearch(atoms, logfile="prerelax.log", trajectory="prerelax.traj")
    dyn.run(fmax=fmax)

    summary = summarize_ase_run(atoms)

    return summary


def loose_relax_positions(
    atoms: Atoms,
    preset: str = "QMOFSet",
    swaps: Dict[str, Any] = None,
) -> Dict[str, Any]:
    """
    Position relaxation with default ENCUT and coarse k-point grid.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    Dict
        Summary of the run.
    """
    swaps = swaps or {}
    defaults = {
        "auto_kpts": {"grid_density": 100},
        "ediff": 1e-4,
        "ediffg": -0.05,
        "encut": None,
        "isif": 2,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 250,
    }
    flags = merge_dicts(defaults, swaps)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    atoms = run_calc(atoms)

    summary = summarize_run(atoms)

    return summary


def loose_relax_volume(
    atoms: Atoms,
    preset: str = "QMOFSet",
    swaps: Dict[str, Any] = None,
) -> Dict[str, Any]:
    """
    Optional: volume relaxation with coarse k-point grid.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    Dict
        Summary of the run.
    """
    swaps = swaps or {}
    defaults = {
        "auto_kpts": {"grid_density": 100},
        "ediffg": -0.03,
        "ibrion": 2,
        "isif": 3,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 500,
    }
    flags = merge_dicts(defaults, swaps)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    atoms = run_calc(atoms)

    summary = summarize_run(atoms)

    return summary


def double_relax(
    atoms: Atoms,
    preset: str = "QMOFSet",
    swaps: Dict[str, Any] = None,
    volume_relax: bool = True,
) -> Dict[str, Any]:
    """
    Double relaxation using production-quality settings.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    volume_relax
        True if a volume relaxation should be performed.

    Returns
    -------
    Dict
        Summary of the run.
    """
    swaps = swaps or {}
    defaults = {
        "ediffg": -0.03,
        "ibrion": 2,
        "isif": 3 if volume_relax else 2,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 500 if volume_relax else 250,
    }

    # Run first relaxation
    flags = merge_dicts(defaults, swaps)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    kpts1 = atoms.calc.kpts
    atoms = run_calc(atoms)

    # Update atoms for
    summary1 = summarize_run(atoms, additional_fields={"name": "relax1"})
    atoms = summary1["atoms"]

    # Reset LREAL
    del defaults["lreal"]

    # Run second relaxation
    flags = merge_dicts(defaults, swaps)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    kpts2 = atoms.calc.kpts

    # Use ISTART = 0 if this goes from vasp_gam --> vasp_std
    if kpts1 == [1, 1, 1] and kpts2 != [1, 1, 1]:
        atoms.calc.set(istart=0)

    atoms = run_calc(atoms)
    summary2 = summarize_run(atoms)

    return {"relax1": summary1, "relax2": summary2}


def static(
    atoms: Atoms,
    preset: str = "QMOFSet",
    swaps: Dict[str, Any] = None,
) -> Dict[str, Any]:
    """
    Static calculation using production-quality settings.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.

    Returns
    -------
    Dict
        Summary of the run.
    """
    swaps = swaps or {}
    defaults = {
        "laechg": True,
        "lcharg": True,
        "lwave": True,
        "nsw": 0,
    }

    # Run static calculation
    flags = merge_dicts(defaults, swaps)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    atoms = run_calc(atoms)

    summary = summarize_run(atoms)

    return summary
