"""Core recipes for GULP"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.gulp import GULP, Conditions
from jobflow import Maker, job

from quacc.schemas.calc import summarize_opt_run, summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_ase_opt, run_calc

GEOM_FILE = GULP.label() + ".got"


@dataclass
class StaticJob(Maker):
    """
    Class to carry out a single-point calculation.
    Note: 'Conditions' are not yet natively supported.

    Parameters
    ----------
    name
        Name of the job.
    gfnff
        True if (periodic) GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    keyword_swaps
        Dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        Dictionary of custom option swap kwargs for the calculator.
    """

    name: str = "GULP-Static"
    gfnff: bool = True
    library: str = None
    keyword_swaps: Dict[str, Any] = field(default_factory=dict)
    option_swaps: Dict[str, Any] = field(default_factory=dict)

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
        default_keywords = {"gfnff": self.gfnff, "gwolf": True if self.gfnff else False}
        default_options = {}

        keywords = merge_dicts(
            default_keywords, self.keyword_swaps, remove_none=True, remove_false=True
        )
        options = merge_dicts(
            default_options, self.option_swaps, remove_none=True, remove_false=True
        )

        gulp_keywords = " ".join(list(keywords.keys()))
        gulp_options = " ".join(list(options.keys()))

        atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options)
        atoms = run_calc(atoms, geom_file=GEOM_FILE)
        summary = summarize_run(
            atoms, input_atoms=atoms, additional_fields={"name": self.name}
        )

        return summary


@dataclass
class RelaxJob(Maker):
    """
    Class to carry out a single-point calculation.
    Note: 'Conditions' are not yet natively supported.

    Parameters
    ----------
    name
        Name of the job.
    gfnff
        True if (periodic) GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    volume_relax
        True if the volume should be relaxed; False if not.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    keyword_swaps
        Dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        Dictionary of custom option swap kwargs for the calculator.
    """

    name: str = "GULP-Static"
    gfnff: bool = True
    library: str = None
    volume_relax: bool = True
    fmax: float = 0.01
    max_steps: int = 1000
    keyword_swaps: Dict[str, Any] = field(default_factory=dict)
    option_swaps: Dict[str, Any] = field(default_factory=dict)

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
        default_keywords = {
            "opti":True,
            "gfnff": self.gfnff,
            "gwolf": True if self.gfnff else False,
            "conp": True if self.volume_relax else False,
            "conv": False if self.volume_relax else True,
        }
        default_options = {}

        keywords = merge_dicts(
            default_keywords, self.keyword_swaps, remove_none=True, remove_false=True
        )
        options = merge_dicts(
            default_options, self.option_swaps, remove_none=True, remove_false=True
        )

        gulp_keywords = " ".join(list(keywords.keys()))
        gulp_options = " ".join(list(options.keys()))

        atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options)
        traj = run_ase_opt(
            atoms,
            fmax=self.fmax,
            max_steps=self.max_steps,
            optimizer="gulp",
            opt_kwargs=self.opt_kwargs,
        )

        summary = summarize_opt_run(
            traj, atoms.calc.parameters, additional_fields={"name": self.name}
        )

        return summary
