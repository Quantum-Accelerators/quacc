"""Core recipes for GULP"""
from __future__ import annotations

import warnings
from typing import Any, Dict

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.gulp import GULP

from quacc.schemas.calc import summarize_run
from quacc.util.basics import merge_dicts
from quacc.util.calc import run_calc


@ct.electron
def StaticJob(
    atoms: Atoms,
    gfnff: bool = True,
    library: str = None,
    keyword_swaps: Dict[str, Any] | None = None,
    option_swaps: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """
    Function to carry out a single-point calculation.
    Note: 'Conditions' are not yet natively supported.

    Parameters
    ----------
    atoms
        .Atoms object
    gfnff
        True if (p)GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    keyword_swaps
        Dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        Dictionary of custom option swap kwargs for the calculator.
    """

    keyword_swaps = keyword_swaps or {}
    option_swaps = option_swaps or {}

    default_keywords = {
        "gfnff": gfnff,
        "gwolf": bool(gfnff and atoms.pbc.any()),
    }
    default_options = {
        "dump every gulp.res": True,
        "output cif gulp.cif": bool(atoms.pbc.any()),
        "output xyz gulp.xyz": not atoms.pbc.any(),
    }

    keywords = merge_dicts(
        default_keywords, keyword_swaps, remove_none=True, remove_false=True
    )
    options = merge_dicts(
        default_options, option_swaps, remove_none=True, remove_false=True
    )

    gulp_keywords = " ".join(list(keywords.keys()))
    gulp_options = list(options.keys())

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options)
    new_atoms = run_calc(atoms, geom_file="gulp.cif" if atoms.pbc.any() else "gulp.xyz")
    summary = summarize_run(new_atoms, input_atoms=atoms)

    return summary


@ct.electron
def RelaxJob(
    atoms: Atoms,
    gfnff: bool = True,
    library: str = None,
    volume_relax: bool = True,
    keyword_swaps: Dict[str, Any] | None = None,
    option_swaps: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """
    Function to carry out a single-point calculation.
    Note: 'Conditions' are not yet natively supported.

    Parameters
    ----------
    atoms
        .Atoms object
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

    keyword_swaps = keyword_swaps or {}
    option_swaps = option_swaps or {}

    if volume_relax and not atoms.pbc.any():
        warnings.warn("Volume relaxation requested but no PBCs found. Ignoring.")
        volume_relax = False

    default_keywords = {
        "opti": True,
        "gfnff": gfnff,
        "gwolf": bool(gfnff and atoms.pbc.any()),
        "conp": bool(volume_relax and atoms.pbc.any()),
        "conv": bool(not volume_relax or not atoms.pbc.any()),
    }
    default_options = {
        "dump every gulp.res": True,
        "output cif gulp.cif": bool(atoms.pbc.any()),
        "output xyz gulp.xyz": not atoms.pbc.any(),
    }

    keywords = merge_dicts(
        default_keywords, keyword_swaps, remove_none=True, remove_false=True
    )
    options = merge_dicts(
        default_options, option_swaps, remove_none=True, remove_false=True
    )

    gulp_keywords = " ".join(list(keywords.keys()))
    gulp_options = list(options.keys())

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options)
    new_atoms = run_calc(atoms, geom_file="gulp.cif" if atoms.pbc.any() else "gulp.xyz")

    if not new_atoms.calc.get_opt_state():
        raise ValueError("Optimization did not converge!")

    summary = summarize_run(new_atoms, input_atoms=atoms)

    return summary
