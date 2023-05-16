"""Core recipes for GULP"""
from __future__ import annotations

import warnings
from copy import deepcopy
from typing import Any

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.gulp import GULP

from quacc.schemas.calc import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import merge_dicts


@ct.electron
def static_job(
    atoms: Atoms,
    gfnff: bool = True,
    library: str | None = None,
    keyword_swaps: dict[str, Any] | None = None,
    option_swaps: dict[str, Any] | None = None,
) -> dict[str, Any]:
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
        dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        dictionary of custom option swap kwargs for the calculator.

    Returns
    -------
    summary
        Dictionary of the calculation summary.
    """

    keyword_swaps = keyword_swaps or {}
    option_swaps = option_swaps or {}
    input_atoms = deepcopy(atoms)

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

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options, library=library)
    atoms = run_calc(atoms, geom_file="gulp.cif" if atoms.pbc.any() else "gulp.xyz")
    summary = summarize_run(
        atoms, input_atoms=input_atoms, additional_fields={"name": "GULP Static"}
    )

    return summary


@ct.electron
def relax_job(
    atoms: Atoms,
    gfnff: bool = True,
    library: str | None = None,
    volume_relax: bool = True,
    keyword_swaps: dict[str, Any] | None = None,
    option_swaps: dict[str, Any] | None = None,
) -> tuple[Atoms, dict[str, Any]]:
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

    Returns
    -------
    summary
        Dictionary of the calculation summary.
    """

    keyword_swaps = keyword_swaps or {}
    option_swaps = option_swaps or {}
    input_atoms = deepcopy(atoms)

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

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options, library=library)
    atoms = run_calc(atoms, geom_file="gulp.cif" if atoms.pbc.any() else "gulp.xyz")

    if not atoms.calc.get_opt_state():
        raise ValueError("Optimization did not converge!")

    summary = summarize_run(
        atoms, input_atoms=input_atoms, additional_fields={"name": "GULP Relax"}
    )

    return summary
