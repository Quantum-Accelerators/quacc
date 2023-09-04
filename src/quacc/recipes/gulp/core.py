"""Core recipes for GULP"""
from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from ase.calculators.gulp import GULP

from quacc import job
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import summarize_run
from quacc.utils.calc import run_calc
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import RunSchema

GEOM_FILE_PBC = "gulp.cif"
GEOM_FILE_NOPBC = "gulp.xyz"


@job
def static_job(
    atoms: Atoms | dict,
    use_gfnff: bool = True,
    library: str | None = None,
    keyword_swaps: dict | None = None,
    option_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation. Note: 'Conditions' are not yet
    natively supported.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    use_gfnff
        True if (p)GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    keyword_swaps
        Dictionary of custom keyword swap kwargs for the calculator. Overrides
        the following defaults:

        ```python
        {
            "gfnff": True if use_gfnff else None,
            "gwolf": True if use_gfnff and atoms.pbc.any() else None,
        }
        ```
    option_swaps
        Dictionary of custom option swap kwargs for the calculator. Overrides
        the following defaults:

        ```python
        {
            "dump every gulp.res": True,
            f"output cif {GEOM_FILE_PBC}": True if atoms.pbc.any() else None,
            f"output xyz {GEOM_FILE_NOPBC}": None if atoms.pbc.any() else True,
        }
        ```

    Returns
    -------
    RunSchema
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """
    atoms = fetch_atoms(atoms)
    keyword_swaps = keyword_swaps or {}
    option_swaps = option_swaps or {}

    default_keywords = {
        "gfnff": True if use_gfnff else None,
        "gwolf": True if use_gfnff and atoms.pbc.any() else None,
    }
    default_options = {
        "dump every gulp.res": True,
        f"output cif {GEOM_FILE_PBC}": True if atoms.pbc.any() else None,
        f"output xyz {GEOM_FILE_NOPBC}": None if atoms.pbc.any() else True,
    }

    keywords = merge_dicts(default_keywords, keyword_swaps)
    options = merge_dicts(default_options, option_swaps)

    gulp_keywords = " ".join(list(keywords.keys()))
    gulp_options = list(options.keys())

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options, library=library)
    final_atoms = run_calc(
        atoms,
        geom_file=GEOM_FILE_PBC if atoms.pbc.any() else GEOM_FILE_NOPBC,
        copy_files=copy_files,
    )

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "GULP Static"},
    )


@job
def relax_job(
    atoms: Atoms | dict,
    use_gfnff: bool = True,
    library: str | None = None,
    relax_cell: bool = False,
    keyword_swaps: dict | None = None,
    option_swaps: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation. Note: 'Conditions' are not yet
    natively supported.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    use_gfnff
        True if (p)GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    relax_cell
        True if the volume should be relaxed; False if not.
    keyword_swaps
        Dictionary of custom keyword swap kwargs for the calculator. Overrides
        the following defaults:

        ```python
        {
            "opti": True,
            "gfnff": True if use_gfnff else None,
            "gwolf": True if use_gfnff and atoms.pbc.any() else None,
            "conp": True if relax_cell and atoms.pbc.any() else None,
            "conv": None if relax_cell and atoms.pbc.any() else True,
        }
        ```
    option_swaps
        Dictionary of custom option swap kwargs for the calculator. Overrides
        the following defaults:

        ```python
        {
            "dump every gulp.res": True,
            f"output cif {GEOM_FILE_PBC}": True if atoms.pbc.any() else None,
            f"output xyz {GEOM_FILE_NOPBC}": None if atoms.pbc.any() else True,
        }
        ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    dict
        Dictionary of results from `quacc.schemas.ase.summarize_run`
    """
    atoms = fetch_atoms(atoms)
    keyword_swaps = keyword_swaps or {}
    option_swaps = option_swaps or {}

    if relax_cell and not atoms.pbc.any():
        warnings.warn(
            "Volume relaxation requested but no PBCs found. Ignoring.", UserWarning
        )
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
        f"output cif {GEOM_FILE_PBC}": True if atoms.pbc.any() else None,
        f"output xyz {GEOM_FILE_NOPBC}": None if atoms.pbc.any() else True,
    }

    keywords = merge_dicts(default_keywords, keyword_swaps)
    options = merge_dicts(default_options, option_swaps)

    gulp_keywords = " ".join(list(keywords.keys()))
    gulp_options = list(options.keys())

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options, library=library)
    final_atoms = run_calc(
        atoms,
        geom_file=GEOM_FILE_PBC if atoms.pbc.any() else GEOM_FILE_NOPBC,
        copy_files=copy_files,
    )

    if not final_atoms.calc.get_opt_state():
        msg = "Optimization did not converge."
        raise ValueError(msg)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "GULP Relax"},
    )
