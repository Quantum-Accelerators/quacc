"""Core recipes for espresso."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.optimize import BFGSLineSearch

from quacc import job
from quacc.atoms.core import check_is_metal
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.recipes.espresso._base import run_and_summarize, run_and_summarize_opt

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from src.quacc.types import EspressoBaseSet

    from quacc.types import Filenames, OptParams, RunSchema, SourceDirectory

BASE_SET_METAL: EspressoBaseSet = {
    "input_data": {
        "system": {"occupations": "smearing", "smearing": "cold", "degauss": 0.01},
        "electrons": {"conv_thr": 1e-8, "mixing_mode": "local-TF", "mixing_beta": 0.35},
        "control": {},
    },
    "kspacing": 0.033,
}

BASE_SET_NON_METAL: EspressoBaseSet = {
    "input_data": {
        "system": {"occupations": "smearing", "smearing": "gaussian", "degauss": 0.005},
        "electrons": {"conv_thr": 1e-8, "mixing_mode": "local-TF", "mixing_beta": 0.35},
    },
    "kspacing": 0.045,
}


@job
def static_job(
    atoms: Atoms,
    preset: str | None = "sssp_1.3.0_pbe_efficiency",
    test_run: bool = False,
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
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
    test_run
        If True, a test run is performed to check that the calculation input_data is correct or
        to generate some files/info if needed.
    copy_files
        Source directory or directories to copy files from. If a `SourceDirectory` or a
        list of `SourceDirectory` is provided, this interface will automatically guess
        which files have to be copied over by looking at the binary and `input_data`.
        If a dict is provided, the mode is manual, keys are source directories and values
        are relative path to files or directories to copy. Glob patterns are supported.
    prev_outdir
        The output directory of a previous calculation. If provided, Quantum Espresso
        will directly read the necessary files from this directory, eliminating the need
        to manually copy files. The directory will be ungzipped if necessary.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    is_metal = check_is_metal(atoms)

    calc_defaults = BASE_SET_METAL if is_metal else BASE_SET_NON_METAL
    calc_defaults["input_data"]["control"] = {"calculation": "scf"}

    return run_and_summarize(
        atoms,
        preset=preset,
        template=EspressoTemplate("pw", test_run=test_run, outdir=prev_outdir),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "pw.x Static"},
        copy_files=copy_files,
    )


@job
def relax_job(
    atoms: Atoms,
    preset: str | None = "sssp_1.3.0_pbe_efficiency",
    relax_cell: bool = False,
    test_run: bool = False,
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
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
    test_run
        If True, a test run is performed to check that the calculation input_data is correct or
        to generate some files/info if needed.
    copy_files
        Source directory or directories to copy files from. If a `SourceDirectory` or a
        list of `SourceDirectory` is provided, this interface will automatically guess
        which files have to be copied over by looking at the binary and `input_data`.
        If a dict is provided, the mode is manual, keys are source directories and values
        are relative path to files or directories to copy. Glob patterns are supported.
    prev_outdir
        The output directory of a previous calculation. If provided, Quantum Espresso
        will directly read the necessary files from this directory, eliminating the need
        to manually copy files. The directory will be ungzipped if necessary.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    is_metal = check_is_metal(atoms)

    calc_defaults = BASE_SET_METAL if is_metal else BASE_SET_NON_METAL
    calc_defaults["input_data"]["control"] = {
        "calculation": "vc-relax" if relax_cell else "relax"
    }

    return run_and_summarize(
        atoms,
        preset=preset,
        template=EspressoTemplate("pw", test_run=test_run, outdir=prev_outdir),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "pw.x Relax"},
        copy_files=copy_files,
    )


@job
def ase_relax_job(
    atoms: Atoms,
    preset: str | None = "sssp_1.3.0_pbe_efficiency",
    autorestart: bool = True,
    relax_cell: bool = False,
    opt_params: OptParams | None = None,
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
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
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.Runner.run_opt][].
    copy_files
        Source directory or directories to copy files from. If a `SourceDirectory` or a
        list of `SourceDirectory` is provided, this interface will automatically guess
        which files have to be copied over by looking at the binary and `input_data`.
        If a dict is provided, the mode is manual, keys are source directories and values
        are relative path to files or directories to copy. Glob patterns are supported.
    prev_outdir
        The output directory of a previous calculation. If provided, Quantum Espresso
        will directly read the necessary files from this directory, eliminating the need
        to manually copy files. The directory will be ungzipped if necessary.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    is_metal = check_is_metal(atoms)

    calc_defaults = BASE_SET_METAL if is_metal else BASE_SET_NON_METAL
    calc_defaults["input_data"]["control"] = {
        "calculation": "scf",
        "tstress": relax_cell,
        "tprnfor": True,
    }

    opt_defaults = {"optimizer": BFGSLineSearch, "relax_cell": relax_cell}

    return run_and_summarize_opt(
        atoms,
        preset=preset,
        template=EspressoTemplate("pw", autorestart=autorestart, outdir=prev_outdir),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        opt_defaults=opt_defaults,
        opt_params=opt_params,
        additional_fields={"name": "pw.x ExternalRelax"},
        copy_files=copy_files,
    )


@job
def post_processing_job(
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
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
        Source directory or directories to copy files from. If a `SourceDirectory` or a
        list of `SourceDirectory` is provided, this interface will automatically guess
        which files have to be copied over by looking at the binary and `input_data`.
        If a dict is provided, the mode is manual, keys are source directories and values
        are relative path to files or directories to copy. Glob patterns are supported.
    prev_outdir
        The output directory of a previous calculation. If provided, Quantum Espresso
        will directly read the necessary files from this directory, eliminating the need
        to manually copy files. The directory will be ungzipped if necessary.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
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

    return run_and_summarize(
        template=EspressoTemplate("pp", test_run=test_run, outdir=prev_outdir),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "pp.x post-processing"},
        copy_files=copy_files,
    )


@job
def non_scf_job(
    atoms: Atoms,
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
    preset: str | None = "sssp_1.3.0_pbe_efficiency",
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
        Source directory or directories to copy files from. If a `SourceDirectory` or a
        list of `SourceDirectory` is provided, this interface will automatically guess
        which files have to be copied over by looking at the binary and `input_data`.
        If a dict is provided, the mode is manual, keys are source directories and values
        are relative path to files or directories to copy. Glob patterns are supported.
    prev_outdir
        The output directory of a previous calculation. If provided, Quantum Espresso
        will directly read the necessary files from this directory, eliminating the need
        to manually copy files. The directory will be ungzipped if necessary.
    preset
        The name of a YAML file containing a list of parameters to use as
        a "preset" for the calculator. quacc will automatically look in the
        `ESPRESSO_PRESET_DIR` (default: quacc/calculators/espresso/presets).
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
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    calc_defaults = {
        "input_data": {"control": {"calculation": "nscf"}},
        "kspacing": 0.033,
    }

    return run_and_summarize(
        atoms,
        preset=preset,
        template=EspressoTemplate("pw", test_run=test_run, outdir=prev_outdir),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "pw.x Non SCF"},
        copy_files=copy_files,
    )
