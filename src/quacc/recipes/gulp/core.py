"""Core recipes for GULP."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.recipes.gulp._base import base_fn

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema


@job
def static_job(
    atoms: Atoms,
    use_gfnff: bool = True,
    keywords: dict[str, Any] | None = None,
    options: dict[str, Any] | None = None,
    library: str | None = None,
) -> RunSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    use_gfnff
        True if (p)GFN-FF should be used; False if not.
    keywords
        Dictionary of custom `keyword` kwargs for the GULP calculator. Set a
        value to `None` to remove a pre-existing key entirely. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    options
        Dictionary of custom `options` kwargs for the GULP calculator. Set a
        value to `None` to remove a pre-existing key entirely. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    library
        Filename of the potential library file, if required.

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

    return base_fn(
        atoms,
        library=library,
        keyword_defaults=keyword_defaults,
        option_defaults=option_defaults,
        keyword_swaps=keywords,
        option_swaps=options,
        additional_fields={"name": "GULP Static"},
    )


@job
def relax_job(
    atoms: Atoms,
    use_gfnff: bool = True,
    relax_cell: bool = False,
    keywords: dict[str, Any] | None = None,
    options: dict[str, Any] | None = None,
    library: str | None = None,
) -> RunSchema:
    """
    Carry out a structure relaxation.

    Parameters
    ----------
    atoms
        Atoms object
    use_gfnff
        True if (p)GFN-FF should be used; False if not.
    relax_cell
        True if the volume should be relaxed; False if not.
    keywords
        Dictionary of custom `keyword` kwargs for the GULP calculator. Set a
        value to `None` to remove a pre-existing key entirely. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    options
        Dictionary of custom `options` kwargs for the GULP calculator. Set a
        value to `None` to remove a pre-existing key entirely. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    library
        Filename of the potential library file, if required.

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

    return base_fn(
        atoms,
        library=library,
        keyword_defaults=keyword_defaults,
        option_defaults=option_defaults,
        keyword_swaps=keywords,
        option_swaps=options,
        additional_fields={"name": "GULP Relax"},
    )
