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
        True if a volume relaxation (ISIF = 3) should be performed.
        False if only the positions (ISIF = 2) should be updated.
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
        defaults = {
            "auto_kpts": {"grid_density": 100},
            "encut": None,
            "ismear": 0,
            "isym": 0,
            "lcharg": False,
            "lreal": "auto",
            "lwave": True,
            "nsw": 0,
            "sigma": 0.05,
        }
        flags = merge_dicts(defaults, swaps, remove_none=True)
        atoms = SmartVasp(atoms, preset=self.preset, **flags)
        dyn = BFGSLineSearch(atoms, logfile="prerelax.log", trajectory="prerelax.traj")
        dyn.run(fmax=5.0)

        # 2. Position relaxation (loose)
        defaults = {
            "auto_kpts": {"grid_density": 100},
            "encut": None,
            "ediff": 1e-5,
            "ediffg": -0.05,
            "isif": 2,
            "ibrion": 2,
            "ismear": 0,
            "isym": 0,
            "lcharg": False,
            "lreal": "auto",
            "lwave": True,
            "nsw": 200,
            "sigma": 0.05,
        }
        flags = merge_dicts(defaults, swaps, remove_none=True)
        atoms = SmartVasp(atoms, preset=self.preset, **flags)
        atoms = run_calc(atoms)

        # 3. Optional: Volume relaxation (loose)
        defaults = {
            "auto_kpts": {"grid_density": 100},
            "ediff": 1e-6,
            "ediffg": -0.02,
            "isif": 2,
            "ibrion": 2,
            "ismear": 0,
            "isym": 0,
            "lcharg": False,
            "lreal": "auto",
            "lwave": True,
            "nsw": 200,
            "sigma": 0.05,
        }
        flags = merge_dicts(defaults, swaps, remove_none=True)
        atoms = SmartVasp(atoms, preset=self.preset, **flags)
        atoms = run_calc(atoms)

        # 4. Double Relaxation
        # This is done for two reasons: a) because it can resolve repadding
        # issues when dV is large; b) because we can use LREAL = Auto for the
        # first relaxation and the default LREAL for the second.
        defaults = {
            "ediff": 1e-6,
            "ediffg": -0.02,
            "isif": 3 if self.volume_relax else 2,
            "ibrion": 2,
            "ismear": 0,
            "isym": 0,
            "lcharg": False,
            "lreal": "auto",
            "lwave": True,
            "nsw": 200,
            "sigma": 0.05,
        }
        flags = merge_dicts(defaults, swaps, remove_none=True)
        old_calc = atoms.calc
        atoms = SmartVasp(atoms, preset=self.preset, **flags)

        # You can't restart a vasp_std calculation from a vasp_gam WAVECAR.
        # Here, we check if we change from vasp_gam to vasp_std and set
        # ISTART = 0 if needed.
        if atoms.calc.kpts != [1, 1, 1] and old_calc.kpts == [1, 1, 1]:
            defaults["istart"] = 0
        atoms = run_calc(atoms)

        # Reset LREAL and ISTART to their default values.
        del defaults["lreal"]
        defaults.pop("istart", None)
        flags = merge_dicts(defaults, swaps, remove_none=True)
        atoms = SmartVasp(atoms, preset=self.preset, **flags)
        atoms = run_calc(atoms)

        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary
