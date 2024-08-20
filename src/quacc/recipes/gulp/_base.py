"""Base jobs for GULP."""

from __future__ import annotations

import os
from logging import getLogger
from typing import TYPE_CHECKING

from ase.calculators.gulp import GULP

from quacc import get_settings
from quacc.runners.ase import Runner
from quacc.schemas.ase import Summarize
from quacc.utils.lists import merge_list_params

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import Filenames, RunSchema, SourceDirectory

LOGGER = getLogger(__name__)

GEOM_FILE_PBC = "gulp.cif"
GEOM_FILE_NOPBC = "gulp.xyz"


def run_and_summarize(
    atoms: Atoms,
    library: str | None = None,
    keyword_defaults: list[str] | None = None,
    option_defaults: list[str] | None = None,
    keyword_swaps: list[str] | None = None,
    option_swaps: list[str] | None = None,
    additional_fields: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
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
    **calc_kwargs
        Any other keyword arguments to pass to the `GULP` calculator.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][]
    """
    keyword_defaults = keyword_defaults or []
    settings = get_settings()
    gulp_cmd = f"{settings.GULP_CMD} < gulp.gin > gulp.got"

    if not atoms.pbc.any():
        if "opti" in keyword_defaults and "conv" not in keyword_defaults:
            keyword_defaults += ["conv"]
        keyword_defaults = [k for k in keyword_defaults if k not in ["gwolf", "conp"]]

    if option_defaults is None:
        option_defaults = []

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

    if settings.GULP_LIB:
        os.environ["GULP_LIB"] = str(settings.GULP_LIB)
    calc = GULP(
        command=gulp_cmd,
        keywords=gulp_keywords,
        options=gulp_options,
        library=library,
        **calc_kwargs,
    )
    final_atoms = Runner(atoms, calc, copy_files=copy_files).run_calc(
        geom_file=GEOM_FILE_PBC if atoms.pbc.any() else GEOM_FILE_NOPBC
    )

    if (
        settings.CHECK_CONVERGENCE
        and "opti" in gulp_keywords
        and not final_atoms.calc.get_opt_state()
    ):
        msg = "Optimization did not converge."
        raise RuntimeError(msg)

    return Summarize(additional_fields=additional_fields).run(final_atoms, atoms)
