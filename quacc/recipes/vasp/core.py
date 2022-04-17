"""Core recipes for VASP"""
from dataclasses import dataclass, field
from typing import Any, Dict

from ase.atoms import Atoms
from jobflow import Maker, job

from quacc.calculators.vasp import Vasp
from quacc.schemas.vasp import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc

GEOM_FILE = "CONTCAR"


@dataclass
class StaticJob(Maker):
    """
    Class to carry out a single-point calculation.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "VASP-Static"
    preset: str = None
    swaps: Dict[str, Any] = field(default_factory=dict)

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
        defaults = {
            "ismear": -5,
            "laechg": True,
            "lcharg": True,
            "lwave": True,
            "nedos": 5001,
            "nsw": 0,
        }
        flags = merge_dicts(defaults, self.swaps)

        calc = Vasp(atoms, preset=self.preset, **flags)
        atoms.calc = calc
        atoms = run_calc(atoms, geom_file=GEOM_FILE)
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary


@dataclass
class RelaxJob(Maker):
    """
    Class to relax a structure.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use.
    volume_relax
        True if a volume relaxation (ISIF = 3) should be performed.
        False if only the positions (ISIF = 2) should be updated.
    swaps
        Dictionary of custom kwargs for the calculator.
    """

    name: str = "VASP-Relax"
    preset: str = None
    volume_relax: bool = True
    swaps: Dict[str, Any] = field(default_factory=dict)

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
        defaults = {
            "ediffg": -0.02,
            "isif": 3 if self.volume_relax else 2,
            "ibrion": 2,
            "isym": 0,
            "lcharg": False,
            "lwave": False,
            "nsw": 200,
        }
        flags = merge_dicts(defaults, self.swaps)

        calc = Vasp(atoms, preset=self.preset, **flags)
        atoms.calc = calc
        atoms = run_calc(atoms, geom_file=GEOM_FILE)
        summary = summarize_run(atoms, additional_fields={"name": self.name})

        return summary


@dataclass
class DoubleRelaxJob(Maker):
    """
    Class to double-relax a structure. This is particularly useful for
    a few reasons:
    1. To carry out a cheaper pre-relaxation before the high-quality run.
    2. To carry out a GGA calculation before a meta-GGA or hybrid calculation
    that requies the GGA wavefunction.
    3. To carry out volume relaxations where large changes in volume
    can require a second relaxation to resolve forces.

    Parameters
    ----------
    name
        Name of the job.
    preset
        Preset to use.
    volume_relax
        True if a volume relaxation (ISIF = 3) should be performed.
        False if only the positions (ISIF = 2) should be updated.
    swaps1
        Dictionary of custom kwargs for the first relaxation.
    swaps2
        Dictionary of custom kwargs for the second relaxation.
    """

    name: str = "VASP-DoubleRelax"
    preset: str = None
    volume_relax: bool = True
    swaps1: Dict[str, Any] = field(default_factory=dict)
    swaps2: Dict[str, Any] = field(default_factory=dict)

    @job
    def make(self, atoms: Atoms) -> Dict[Dict[str, Any], Dict[str, Any]]:
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
            "ediffg": -0.02,
            "isif": 3 if self.volume_relax else 2,
            "ibrion": 2,
            "isym": 0,
            "lcharg": False,
            "lwave": True,
            "nsw": 200,
        }

        # Run first relaxation
        flags = merge_dicts(defaults, self.swaps1)
        calc = Vasp(atoms, preset=self.preset, **flags)
        atoms.calc = calc
        kpts1 = atoms.calc.kpts
        atoms = run_calc(atoms, geom_file=GEOM_FILE)
        summary1 = summarize_run(atoms, additional_fields={"name": self.name})

        # Run second relaxation
        flags = merge_dicts(defaults, self.swaps2)
        calc = Vasp(summary1["atoms"], preset=self.preset, **flags)
        atoms.calc = calc
        kpts2 = atoms.calc.kpts

        # Use ISTART = 0 if this goes from vasp_gam --> vasp_std
        if kpts1 == [1, 1, 1] and kpts2 != [1, 1, 1]:
            atoms.calc.set(istart=0)

        atoms = run_calc(atoms, geom_file=GEOM_FILE)
        summary2 = summarize_run(atoms, additional_fields={"name": self.name})

        return {"relax1": summary1, "relax2": summary2}
