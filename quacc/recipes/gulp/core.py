"""Core recipes for GULP"""
from __future__ import annotations

import warnings
from copy import deepcopy

import covalent as ct
from ase.atoms import Atoms
from ase.calculators.gulp import GULP

from quacc.schemas.ase import summarize_run
from quacc.util.calc import run_calc
from quacc.util.dicts import remove_dict_empties


@ct.electron
def static_job(
    atoms: Atoms,
    use_gfnff: bool = True,
    library: str | None = None,
    keyword_swaps: dict | None = None,
    option_swaps: dict | None = None,
) -> dict:
    """
    Carry out a single-point calculation.
    Note: 'Conditions' are not yet natively supported.

    Parameters
    ----------
    atoms
        Atoms object
    use_gfnff
        True if (p)GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    keyword_swaps
        dictionary of custom keyword swap kwargs for the calculator.
            default_keywords = {
                "gfnff": True if use_gfnff else None,
                "gwolf": True if use_gfnff and atoms.pbc.any() else None,
            }
    option_swaps
        dictionary of custom option swap kwargs for the calculator.
            default_options = {
                "dump every gulp.res": True,
                "output cif gulp.cif": True if atoms.pbc.any() else None,
                "output xyz gulp.xyz": None if atoms.pbc.any() else True,
            }

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """

    keyword_swaps = keyword_swaps or {}
    option_swaps = option_swaps or {}
    input_atoms = deepcopy(atoms)

    default_keywords = {
        "gfnff": True if use_gfnff else None,
        "gwolf": True if use_gfnff and atoms.pbc.any() else None,
    }
    default_options = {
        "dump every gulp.res": True,
        "output cif gulp.cif": True if atoms.pbc.any() else None,
        "output xyz gulp.xyz": None if atoms.pbc.any() else True,
    }

    keywords = remove_dict_empties(default_keywords | keyword_swaps)
    options = remove_dict_empties(default_options | option_swaps)

    gulp_keywords = " ".join(list(keywords.keys()))
    gulp_options = list(options.keys())

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options, library=library)
    atoms = run_calc(atoms, geom_file="gulp.cif" if atoms.pbc.any() else "gulp.xyz")

    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        additional_fields={"name": "GULP Static"},
    )


@ct.electron
def relax_job(
    atoms: Atoms,
    use_gfnff: bool = True,
    library: str | None = None,
    relax_cell: bool = True,
    keyword_swaps: dict | None = None,
    option_swaps: dict | None = None,
) -> dict:
    """
    Carry out a single-point calculation.
    Note: 'Conditions' are not yet natively supported.

    Parameters
    ----------
    atoms
        Atoms object
    use_gfnff
        True if (p)GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    relax_cell
        True if the volume should be relaxed; False if not.
    keyword_swaps
        Dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        Dictionary of custom option swap kwargs for the calculator.

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """

    keyword_swaps = keyword_swaps or {}
    option_swaps = option_swaps or {}
    input_atoms = deepcopy(atoms)

    if relax_cell and not atoms.pbc.any():
        warnings.warn("Volume relaxation requested but no PBCs found. Ignoring.")
        relax_cell = False

    default_keywords = {
        "opti": True,
        "gfnff": True if use_gfnff else None,
        "gwolf": True if use_gfnff and atoms.pbc.any() else None,
        "conp": True if relax_cell and atoms.pbc.any() else None,
        "conv": None if relax_cell and atoms.pbc.any() else True,
    }
    default_options = {
        "dump every gulp.res": True,
        "output cif gulp.cif": True if atoms.pbc.any() else None,
        "output xyz gulp.xyz": None if atoms.pbc.any() else True,
    }

    keywords = remove_dict_empties(default_keywords | keyword_swaps)
    options = remove_dict_empties(default_options | option_swaps)

    gulp_keywords = " ".join(list(keywords.keys()))
    gulp_options = list(options.keys())

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options, library=library)
    atoms = run_calc(atoms, geom_file="gulp.cif" if atoms.pbc.any() else "gulp.xyz")

    if not atoms.calc.get_opt_state():
        raise ValueError("Optimization did not converge!")

    return summarize_run(
        atoms,
        input_atoms=input_atoms,
        additional_fields={"name": "GULP Relax"},
    )
