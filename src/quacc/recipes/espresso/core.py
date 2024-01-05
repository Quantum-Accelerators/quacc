"""Core recipes for espresso."""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.recipes.espresso._base import base_fn

if TYPE_CHECKING:
    from pathlib import Path

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema


@job
def static_job(
    atoms: Atoms,
    preset: str | None = "sssp_1.3.0_pbe_efficiency",
    parallel_info: dict[str] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
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
        `ESPRESSO_PRESET_DIR` (default: quacc/calculators/espresso/presets).
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
        Additional keyword arguments to pass to the Espresso calculator. See the
        docstring of [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = {"input_data": {"control": {"calculation": "scf"}}}

    return base_fn(
        atoms,
        preset=preset,
        template=EspressoTemplate("pw"),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "pw.x Static"},
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    preset: str | None = "sssp_1.3.0_pbe_efficiency",
    relax_cell: bool = False,
    parallel_info: dict[str] | None = None,
    copy_files: str | Path | list[str | Path] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a structure relaxation with pw.x.

    Parameters
    ----------
    atoms
        The Atoms object.
    preset
        The name of a YAML file containing a list of parameters to use as
        a "preset" for the calculator. quacc will automatically look in the
        `ESPRESSO_PRESET_DIR` (default: quacc/calculators/espresso/presets).
    relax_cell
        Whether to relax the cell or not.
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
        Additional keyword arguments to pass to the Espresso calculator. See the
        docstring of [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = {
        "input_data": {
            "control": {"calculation": "vc-relax" if relax_cell else "relax"}
        }
    }

    return base_fn(
        atoms,
        preset=preset,
        template=EspressoTemplate("pw"),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "pw.x Relax"},
        copy_files=copy_files,
    )


@job
def post_processing_job(
    prev_dir: str | Path, parallel_info: dict[str] | None = None, **calc_kwargs
) -> RunSchema:
    """
    Function to carry out a basic pp.x calculation (post-processing).
    It is mainly used to extract the charge density from a previous pw.x calculation.
    and perform simple to complex post-processing on it. Fore more details please see
    (https://www.quantum-espresso.org/Doc/INPUT_PP.html)

    Parameters
    ----------
    prev_dir
        Outdir of the previously ran pw.x calculation. This is used to copy
        the entire tree structure of that directory to the working directory
        of this calculation.
    parallel_info
        Dictionary containing information about the parallelization of the
        calculation. See the ASE documentation for more information.
    **calc_kwargs
        calc_kwargs dictionary possibly containing the following keys:

        - input_data: dict
        - additional_fields: list[str] | str

        See the docstring of ase.io.espresso.write_fortran_namelist for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    calc_defaults = {
        "input_data": {
            "inputpp": {"plot_num": 0},
            "plot": {
                "iflag": 3,
                "output_format": 6,
                "fileout": "pseudo_charge_density.cube",
            },
        }
    }

    return base_fn(
        template=EspressoTemplate("pp"),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "pp.x post-processing"},
        copy_files=prev_dir,
    )
