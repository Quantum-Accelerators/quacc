"""Core recipes for GULP"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from ase.calculators.gulp import GULP

from quacc import SETTINGS, job
from quacc.runners.calc import run_calc
from quacc.schemas import fetch_atoms
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from ase import Atoms

    from quacc.schemas.ase import RunSchema

logger = logging.getLogger(__name__)

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
    Carry out a single-point calculation.

    ??? Note

        Keyword Defaults:

        ```python
        {
            "gfnff": True if use_gfnff else None,
            "gwolf": True if use_gfnff and atoms.pbc.any() else None,
        }
        ```

        Option Defaultss:

        ```python
        {
            "dump every gulp.res": True,
            f"output cif {GEOM_FILE_PBC}": True if atoms.pbc.any() else None,
            f"output xyz {GEOM_FILE_NOPBC}": None if atoms.pbc.any() else True,
        }
        ```

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
        Dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        Dictionary of custom option swap kwargs for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    keyword_defaults = {
        "gfnff": True if use_gfnff else None,
        "gwolf": True if use_gfnff else None,
    }
    option_defaults = {"dump every gulp.res": True}

    return _base_job(
        atoms,
        library=library,
        keyword_defaults=keyword_defaults,
        option_defaults=option_defaults,
        keyword_swaps=keyword_swaps,
        option_swaps=option_swaps,
        additional_fields={"name": "GULP Static"},
        copy_files=copy_files,
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
    Carry out a structure relaxation.

    ??? Note

        Keyword Defaults:

        ```python
        {
            "opti": True,
            "gfnff": True if use_gfnff else None,
            "gwolf": True if use_gfnff and atoms.pbc.any() else None,
            "conp": True if relax_cell and atoms.pbc.any() else None,
            "conv": None if relax_cell and atoms.pbc.any() else True,
        }
        ```

        Option Defaults:

        ```python
        {
            "dump every gulp.res": True,
            f"output cif {GEOM_FILE_PBC}": True if atoms.pbc.any() else None,
            f"output xyz {GEOM_FILE_NOPBC}": None if atoms.pbc.any() else True,
        }
        ```

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
        Dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        Dictionary of custom option swap kwargs for the calculator.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    dict
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    keyword_defaults = {
        "opti": True,
        "gfnff": True if use_gfnff else None,
        "gwolf": True if use_gfnff else None,
        "conp": True if relax_cell else None,
        "conv": None if relax_cell else True,
    }
    option_defaults = {"dump every gulp.res": True}

    return _base_job(
        atoms,
        library=library,
        keyword_defaults=keyword_defaults,
        option_defaults=option_defaults,
        keyword_swaps=keyword_swaps,
        option_swaps=option_swaps,
        additional_fields={"name": "GULP Relax"},
        copy_files=copy_files,
    )


def _base_job(
    atoms: Atoms | dict,
    library: str | None = None,
    keyword_defaults: dict | None = None,
    option_defaults: dict | None = None,
    keyword_swaps: dict | None = None,
    option_swaps: dict | None = None,
    additional_fields: dict | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Base job function for GULP recipes.

    Parameters
    ----------
    atoms
        Atoms object or a dictionary with the key "atoms" and an Atoms object as
        the value
    library
        Filename of the potential library file, if required.
    keyword_defaults
        Default keywords for calculator.
    option_defaults
        Default options for calculator.
    keyword_swaps
        Dictionary of custom keyword swap kwargs for the calculator.
    option_swaps
        Dictionary of custom option swap kwargs for the calculator.
    additional_fields
        Additional field to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    atoms = fetch_atoms(atoms)
    keyword_defaults = keyword_defaults or {}

    if not atoms.pbc.any():
        if keyword_defaults.get("opti") and not keyword_defaults.get("conv"):
            keyword_defaults["conv"] = True
        for k in ["gwolf", "conp"]:
            keyword_defaults.pop(k, None)

    option_defaults = merge_dicts(
        option_defaults,
        {
            f"output cif {GEOM_FILE_PBC}": True if atoms.pbc.any() else None,
            f"output xyz {GEOM_FILE_NOPBC}": None if atoms.pbc.any() else True,
        },
    )

    keywords = merge_dicts(keyword_defaults, keyword_swaps)
    options = merge_dicts(option_defaults, option_swaps)

    gulp_keywords = " ".join(list(keywords.keys()))
    gulp_options = list(options.keys())

    atoms.calc = GULP(keywords=gulp_keywords, options=gulp_options, library=library)
    final_atoms = run_calc(
        atoms,
        geom_file=GEOM_FILE_PBC if atoms.pbc.any() else GEOM_FILE_NOPBC,
        copy_files=copy_files,
    )

    if (
        SETTINGS.CHECK_CONVERGENCE
        and "opti" in gulp_keywords
        and not final_atoms.calc.get_opt_state()
    ):
        msg = "Optimization did not converge."
        raise ValueError(msg)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields=additional_fields,
    )
