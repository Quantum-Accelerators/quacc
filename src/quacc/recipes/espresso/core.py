"""Core recipes for espresso."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import LBFGS

from quacc import job
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.recipes.espresso._base import base_fn, base_opt_fn

if TYPE_CHECKING:
    from typing import Any

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory


@job
def static_job(
    atoms: Atoms,
    preset: str | None = "sssp_1.3.0_pbe_efficiency",
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
    test_run
        If True, a test run is performed to check that the calculation input_data is correct or
        to generate some files/info if needed.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {"input_data": {"control": {"calculation": "scf"}}}

    return base_fn(
        atoms,
        preset=preset,
        template=EspressoTemplate("pw", test_run=test_run),
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
    test_run: bool = False,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
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
    test_run
        If True, a test run is performed to check that the calculation input_data is correct or
        to generate some files/info if needed.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "input_data": {
            "control": {"calculation": "vc-relax" if relax_cell else "relax"}
        }
    }

    return base_fn(
        atoms,
        preset=preset,
        template=EspressoTemplate("pw", test_run=test_run),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "pw.x Relax"},
        copy_files=copy_files,
    )


@job
def ase_relax_job(
    atoms: Atoms,
    preset: str | None = "sssp_1.3.0_pbe_efficiency",
    autorestart: bool = True,
    relax_cell: bool = False,
    parallel_info: dict[str] | None = None,
    opt_params: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a structure relaxation with pw.x using ASE
    external optimizers.

    Parameters
    ----------
    atoms
        The Atoms object.
    preset
        The name of a YAML file containing a list of parameters to use as
        a "preset" for the calculator. quacc will automatically look in the
        `ESPRESSO_PRESET_DIR` (default: quacc/calculators/espresso/presets).
    autorestart
        Whether to automatically turn on the restart flag after the first
        calculation. This avoids recomputing everything from scratch at each
        step of the optimization.
    relax_cell
        Whether to relax the cell or not.
    parallel_info
        Dictionary containing information about the parallelization of the
        calculation. See the ASE documentation for more information.
    opt_params
        Dictionary of parameters to pass to the optimizer. pass "optimizer"
        to change the optimizer being used. "fmax" and "max_steps" are commonly
        used keywords. See the ASE documentation for more information.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {
        "input_data": {
            "control": {"calculation": "scf", "tstress": relax_cell, "tprnfor": True}
        }
    }

    opt_defaults = {"optimizer": LBFGS}

    return base_opt_fn(
        atoms,
        preset=preset,
        relax_cell=relax_cell,
        template=EspressoTemplate("pw", autorestart=autorestart),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        parallel_info=parallel_info,
        additional_fields={"name": "pw.x ExternalRelax"},
        copy_files=copy_files,
    )


@job
def post_processing_job(
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames],
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic pp.x calculation (post-processing).
    It is mainly used to extract the charge density from a previous pw.x calculation.
    and perform simple to complex post-processing on it. Fore more details please see
    https://www.quantum-espresso.org/Doc/INPUT_PP.html

    Parameters
    ----------
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    parallel_info
        Dictionary containing information about the parallelization of the
        calculation. See the ASE documentation for more information.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
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
        template=EspressoTemplate("pp", test_run=test_run),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "pp.x post-processing"},
        copy_files=copy_files,
    )


@job
def non_scf_job(
    atoms: Atoms,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames],
    preset: str | None = "sssp_1.3.0_pbe_efficiency",
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic NSCF calculation with pw.x.

    Parameters
    ----------
    atoms
        The Atoms object.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    preset
        The name of a YAML file containing a list of parameters to use as
        a "preset" for the calculator. quacc will automatically look in the
        `ESPRESSO_PRESET_DIR` (default: quacc/calculators/espresso/presets).
    parallel_info
        Dictionary containing information about the parallelization of the
        calculation. See the ASE documentation for more information.
    test_run
        If True, a test run is performed to check that the calculation input_data is correct or
        to generate some files/info if needed.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    calc_defaults = {"input_data": {"control": {"calculation": "nscf"}}}

    return base_fn(
        atoms,
        preset=preset,
        template=EspressoTemplate("pw", test_run=test_run),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "pw.x Non SCF"},
        copy_files=copy_files,
    )
