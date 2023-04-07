"""Core recipes for GULP"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any, Dict

from ase.atoms import Atoms
from ase.calculators.gulp import GULP
from jobflow import Maker, job

from quacc.schemas.calc import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc


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
        True if (p)GFN-FF should be used; False if not.
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

        default_keywords = {
            "gfnff": self.gfnff,
            "gwolf": bool(self.gfnff and atoms.pbc.any()),
        }
        default_options = {
            "dump every gulp.res": True,
            "output cif gulp.cif": bool(atoms.pbc.any()),
            "output xyz gulp.xyz": not atoms.pbc.any(),
        }

        keywords = merge_dicts(
            default_keywords, self.keyword_swaps, remove_none=True, remove_false=True
        )
        options = merge_dicts(
            default_options, self.option_swaps, remove_none=True, remove_false=True
        )

        gulp_keywords = " ".join(list(keywords.keys()))
        gulp_options = list(options.keys())

        atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options)
        new_atoms = run_calc(
            atoms, geom_file="gulp.cif" if atoms.pbc.any() else "gulp.xyz"
        )
        summary = summarize_run(
            new_atoms, input_atoms=atoms, additional_fields={"name": self.name}
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
        True if (p)GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    volume_relax
        True if the volume should be relaxed; False if not.
    keyword_swaps
        Dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        Dictionary of custom option swap kwargs for the calculator.
    """

    name: str = "GULP-Relax"
    gfnff: bool = True
    library: str = None
    volume_relax: bool = True
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
        if self.volume_relax and not atoms.pbc.any():
            warnings.warn("Volume relaxation requested but no PBCs found. Ignoring.")
            self.volume_relax = False

        default_keywords = {
            "opti": True,
            "gfnff": self.gfnff,
            "gwolf": bool(self.gfnff and atoms.pbc.any()),
            "conp": bool(self.volume_relax and atoms.pbc.any()),
            "conv": bool(not self.volume_relax or not atoms.pbc.any()),
        }
        default_options = {
            "dump every gulp.res": True,
            "output cif gulp.cif": bool(atoms.pbc.any()),
            "output xyz gulp.xyz": not atoms.pbc.any(),
        }

        keywords = merge_dicts(
            default_keywords, self.keyword_swaps, remove_none=True, remove_false=True
        )
        options = merge_dicts(
            default_options, self.option_swaps, remove_none=True, remove_false=True
        )

        gulp_keywords = " ".join(list(keywords.keys()))
        gulp_options = list(options.keys())

        atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options)
        new_atoms = run_calc(
            atoms, geom_file="gulp.cif" if atoms.pbc.any() else "gulp.xyz"
        )

        if not new_atoms.calc.get_opt_state():
            raise ValueError("Optimization did not converge!")

        summary = summarize_run(
            new_atoms, input_atoms=atoms, additional_fields={"name": self.name}
        )

        return summary
