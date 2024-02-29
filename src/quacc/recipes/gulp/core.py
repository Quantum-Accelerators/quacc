"""Core recipes for GULP."""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.recipes.gulp._base import base_fn

if TYPE_CHECKING:
    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    use_gfnff: bool = True,
    keywords: list[str] | None = None,
    options: list[str] | None = None,
    library: str | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        List of custom `keyword` kwargs for the GULP calculator. To remove entries
        from the defaults, put a `#` in front of the name. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    options
        List of custom `options` kwargs for the GULP calculator. To remove entries
        from the defaults, put a `#` in front of the name. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    library
        Filename of the potential library file, if required.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    keyword_defaults = ["gfnff", "gwolf"] if use_gfnff else []
    option_defaults = ["dump every gulp.res"]

    return base_fn(
        atoms,
        library=library,
        keyword_defaults=keyword_defaults,
        option_defaults=option_defaults,
        keyword_swaps=keywords,
        option_swaps=options,
        additional_fields={"name": "GULP Static"},
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    use_gfnff: bool = True,
    relax_cell: bool = False,
    keywords: list[str] | None = None,
    options: list[str] | None = None,
    library: str | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
        List of custom `keyword` kwargs for the GULP calculator. To remove entries
        from the defaults, put a `#` in front of the name. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    options
        Dictionary of custom `options` kwargs for the GULP calculator. To remove entries
        from the defaults, put a `#` in front of the name. For a list of
        available keys, refer to the `ase.calculators.gulp.GULP` calculator.
    library
        Filename of the potential library file, if required.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    keyword_defaults = ["opti", "conp" if relax_cell else "conv"]
    if use_gfnff:
        keyword_defaults += ["gfnff", "gwolf"]

    option_defaults = ["dump every gulp.res"]

    return base_fn(
        atoms,
        library=library,
        keyword_defaults=keyword_defaults,
        option_defaults=option_defaults,
        keyword_swaps=keywords,
        option_swaps=options,
        additional_fields={"name": "GULP Relax"},
        copy_files=copy_files,
    )
