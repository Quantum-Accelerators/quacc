"""Base jobs for GULP."""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING

from ase.calculators.gulp import GULP

from quacc import SETTINGS
from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.lists import merge_list_params

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory
logger = logging.getLogger(__name__)

GEOM_FILE_PBC = "gulp.cif"
GEOM_FILE_NOPBC = "gulp.xyz"
GULP_CMD = f"{SETTINGS.GULP_CMD} < gulp.gin > gulp.got"


def base_fn(
    atoms: Atoms,
    library: str | None = None,
    keyword_defaults: list[str] | None = None,
    option_defaults: list[str] | None = None,
    keyword_swaps: list[str] | None = None,
    option_swaps: list[str] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        List of custom `keyword` kwargs for the GULP calculator. To remove entries
        from the defaults, put a `#` in front of the name. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    option_swaps
        Dictionary of custom `options` kwargs for the GULP calculator. To remove entries
        from the defaults, put a `#` in front of the name. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    additional_fields
        Additional field to supply to the summarizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """
    keyword_defaults = keyword_defaults or []

    if not atoms.pbc.any():
        if "opti" in keyword_defaults and "conv" not in keyword_defaults:
            keyword_defaults += ["conv"]
        keyword_defaults = [k for k in keyword_defaults if k not in ["gwolf", "conp"]]

    option_defaults += [
        (
            f"output cif {GEOM_FILE_PBC}"
            if atoms.pbc.any()
            else f"output xyz {GEOM_FILE_NOPBC}"
        )
    ]

    keywords = merge_list_params(keyword_defaults, keyword_swaps)
    options = merge_list_params(option_defaults, option_swaps)

    gulp_keywords = " ".join(keywords)
    gulp_options = list(options)

    if SETTINGS.GULP_LIB:
        os.environ["GULP_LIB"] = str(SETTINGS.GULP_LIB)
    atoms.calc = GULP(
        command=GULP_CMD, keywords=gulp_keywords, options=gulp_options, library=library
    )
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
        raise RuntimeError(msg)

    return summarize_run(final_atoms, atoms, additional_fields=additional_fields)
