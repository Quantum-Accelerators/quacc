from dataclasses import dataclass
from typing import Any, Dict

from ase.atoms import Atoms
from ase.optimize import BFGSLineSearch
from jobflow import Maker, job

from quacc.calculators.vasp import SmartVasp
from quacc.schemas.vasp import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc


@dataclass
class MultiRelaxMaker(Maker):
    """
    Class to relax a structure in a multi-step process for increased
    computational efficiency. This is all done in a single compute job.

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

    name: str = "VASP-MultiRelax"
    preset: str = None
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
        atoms = prerelax(atoms, self.preset, swaps, fmax=5.0)

        # 2. Position relaxation (loose)
        atoms = loose_relax_positions(atoms, self.preset, swaps)

        # 3. Optional: Volume relaxation (loose)
        if self.volume_relax:
            atoms = loose_relax_volume(atoms, self.preset, swaps)

        # 4. Double Relaxation
        # This is done for two reasons: a) because it can resolve repadding
        # issues when dV is large; b) because we can use LREAL = Auto for the
        # first relaxation and the default LREAL for the second.
        atoms = double_relax(atoms, self.preset, swaps, volume_relax=self.volume_relax)

        # Make summary of run
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary


def prerelax(
    atoms: Atoms,
    preset: str,
    swaps: Dict[str, Any],
    fmax: float = 5.0,
) -> Atoms:
    """
    A "pre-relaxation" with BFGSLineSearch to resolve very high forces.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use. Applies for all jobs.
    swaps
        Dictionary of custom kwargs for the calculator. Applies for all jobs.
    fmax
        Maximum force in eV/A.

    Returns
    -------
    .Atoms object
    """
    defaults = {
        "auto_kpts": {"grid_density": 100},
        "ediff": 1e-4,
        "encut": None,
        "ismear": 0,
        "isym": 0,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nelm": 225,
        "nsw": 0,
        "sigma": 0.05,
    }
    flags = merge_dicts(defaults, swaps, remove_none=True)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    dyn = BFGSLineSearch(atoms, logfile="prerelax.log", trajectory="prerelax.traj")
    dyn.run(fmax=fmax)

    return atoms


def loose_relax_positions(
    atoms: Atoms,
    preset: str,
    swaps: Dict[str, Any],
) -> Atoms:
    """
    Position relaxation with default ENCUT and coarse k-point grid.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use. Applies for all jobs.
    swaps
        Dictionary of custom kwargs for the calculator. Applies for all jobs.

    Returns
    -------
    .Atoms object
    """
    defaults = {
        "auto_kpts": {"grid_density": 100},
        "encut": None,
        "ediff": 1e-4,
        "ediffg": -0.05,
        "isif": 2,
        "ibrion": 2,
        "ismear": 0,
        "isym": 0,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 250,
        "sigma": 0.05,
    }
    flags = merge_dicts(defaults, swaps, remove_none=True)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    atoms = run_calc(atoms)

    return atoms


def loose_relax_volume(
    atoms: Atoms,
    preset: str,
    swaps: Dict[str, Any],
) -> Atoms:
    """
    Optional: volume relaxation with coarse k-point grid.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use. Applies for all jobs.
    swaps
        Dictionary of custom kwargs for the calculator. Applies for all jobs.

    Returns
    -------
    .Atoms object
    """
    defaults = {
        "auto_kpts": {"grid_density": 100},
        "ediff": 1e-6,
        "ediffg": -0.02,
        "isif": 3,
        "ibrion": 2,
        "ismear": 0,
        "isym": 0,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 500,
        "sigma": 0.05,
    }
    flags = merge_dicts(defaults, swaps, remove_none=True)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    atoms = run_calc(atoms)

    return atoms


def double_relax(
    atoms: Atoms,
    preset: str,
    swaps: Dict[str, Any],
    volume_relax: bool = True,
) -> Atoms:
    """
    Double relaxation using production-quality settings.

    Parameters
    ----------
    atoms
        .Atoms object
    preset
        Preset to use. Applies for all jobs.
    swaps
        Dictionary of custom kwargs for the calculator. Applies for all jobs.
    volume_relax
        True if a volume relaxation should be performed.

    Returns
    -------
    .Atoms object
    """

    defaults = {
        "ediff": 1e-6,
        "ediffg": -0.02,
        "isif": 3 if volume_relax else 2,
        "ibrion": 2,
        "ismear": 0,
        "isym": 0,
        "lcharg": False,
        "lreal": "auto",
        "lwave": True,
        "nsw": 250,
        "sigma": 0.05,
    }
    flags = merge_dicts(defaults, swaps, remove_none=True)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    kpts1 = atoms.calc.kpts
    atoms = run_calc(atoms)

    # Reset LREAL and ISTART to their default values.
    del defaults["lreal"]

    flags = merge_dicts(defaults, swaps, remove_none=True)
    atoms = SmartVasp(atoms, preset=preset, **flags)
    kpts2 = atoms.calc.kpts

    # Use ISTART = 0 if this goes from vasp_gam --> vasp_std
    if kpts1 == [1, 1, 1] and kpts2 != [1, 1, 1]:
        atoms.calc.set(istart=0)

    atoms = run_calc(atoms)

    return atoms
