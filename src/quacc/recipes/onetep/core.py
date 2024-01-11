"""Core recipes for Onetep."""
from __future__ import annotations

from typing import TYPE_CHECKING

from ase import Atoms

from quacc import job
from quacc.recipes.onetep._base import base_fn

if TYPE_CHECKING:
    from typing import Any

    from quacc.schemas._aliases.ase import RunSchema


@job
def static_job(
    atoms: Atoms,
    preset: str | None = None,
    copy_files: list[str] | None = None,
    parallel_info: dict[str] | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic SCF calculation with pw.x.

    Parameters
    ----------
    atoms
        The Atoms object.
    preset
        The name of a YAML file containing a list of parameters to use as
        a "preset" for the calculator. quacc will automatically look in the
        `ONETEP_PRESET_DIR` (default: quacc/calculators/onetep/presets).
    parallel_info
        Dictionary containing information about the parallelization of the
        calculation. See the ASE documentation for more information.
    copy_files
        List of files to copy to the calculation directory. Useful for copying
        files from a previous calculation. This parameter can either be a string
        or a list of strings.

        If a string is provided, it is assumed to be a path to a directory,
        all of the child tree structure of that directory is going to be copied to the
        scratch of this calculation. For phonon_job this is what most users will want to do.

        If a list of strings is provided, each string point to a specific file. In this case
        it is important to note that no directory structure is going to be copied, everything
        is copied at the root of the temporary directory.
    **calc_kwargs
        Additional keyword arguments to pass to the ONETEP calculator. See the
        docstring of [quacc.calculators.onetep.Onetep][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = {
        "keywords": {
            "output_detail": "verbose",
            "do_properties": True,
            "cutoff_energy": "600 eV",
            "task": "SinglePoint",
        }
    }

    return base_fn(
        atoms,
        preset=preset,
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields=additional_fields,
        copy_files=copy_files,
    )
