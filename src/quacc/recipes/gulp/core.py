"""Core recipes for GULP."""
from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING

from ase.calculators.gulp import GULP

from quacc import SETTINGS, job
from quacc.runners.calc import run_ase_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import merge_dicts

if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms

    from quacc.schemas.ase import RunSchema

logger = logging.getLogger(__name__)

GEOM_FILE_PBC = "gulp.cif"
GEOM_FILE_NOPBC = "gulp.xyz"
GULP_CMD = f"{SETTINGS.GULP_CMD} < gulp.gin > gulp.got"


@job
def static_job(
    atoms: Atoms,
    use_gfnff: bool = True,
    library: str | None = None,
    keyword_swaps: dict[str, Any] | None = None,
    option_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    use_gfnff
        True if (p)GFN-FF should be used; False if not.
    library
        Filename of the potential library file, if required.
    keyword_swaps
        Dictionary of custom `keyword` kwargs for the GULP calculator. Set a
        value to `None` to remove a pre-existing key entirely. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.

        !!! Info "Keyword defaults"

            ```python
            {
                "gfnff": True if use_gfnff else None,
                "gwolf": True if use_gfnff and atoms.pbc.any() else None,
            }
            ```
    option_swaps
        Dictionary of custom `options` kwargs for the GULP calculator. Set a
        value to `None` to remove a pre-existing key entirely. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.

        !!! Info "Option defaults"

            ```python
            {"dump every gulp.res": True}
            ```
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
    atoms: Atoms,
    use_gfnff: bool = True,
    library: str | None = None,
    relax_cell: bool = False,
    keyword_swaps: dict[str, Any] | None = None,
    option_swaps: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Carry out a structure relaxation.

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
        Dictionary of custom `keyword` kwargs for the GULP calculator. Set a
        value to `None` to remove a pre-existing key entirely. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.

        !!! Info "Keyword defaults"

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
        Dictionary of custom `options` kwargs for the GULP calculator. Set a
        value to `None` to remove a pre-existing key entirely. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.

        !!! Info "Option defaults"

            ```python
            {"dump every gulp.res": True}
            ```
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
    atoms: Atoms,
    library: str | None = None,
    keyword_defaults: dict[str, Any] | None = None,
    option_defaults: dict[str, Any] | None = None,
    keyword_swaps: dict[str, Any] | None = None,
    option_swaps: dict[str, Any] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: list[str] | None = None,
) -> RunSchema:
    """
    Base job function for GULP recipes.

    Parameters
    ----------
    atoms
        Atoms object
    library
        Filename of the potential library file, if required.
    keyword_defaults
        Default `keywords` for calculator.
    option_defaults
        Default `options` for calculator.
    keyword_swaps
        Dictionary of custom `keyword` kwargs for the GULP calculator. Set a
        value to `None` to remove a pre-existing key entirely. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    option_swaps
        Dictionary of custom `options` kwargs for the GULP calculator. Set a
        value to `None` to remove a pre-existing key entirely. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    additional_fields
        Additional field to supply to the summarizer.
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
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

    if SETTINGS.GULP_LIB:
        os.environ["GULP_LIB"] = str(SETTINGS.GULP_LIB)
    atoms.calc = GULP(
        command=GULP_CMD,
        keywords=gulp_keywords,
        options=gulp_options,
        library=library,
    )
    final_atoms = run_ase_calc(
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
        raise RuntimeError(msg)

    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields=additional_fields,
    )
