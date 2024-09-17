"""
This module, 'bands.py', contains recipes for performing bands and fermi surface calculations using the
bands.x and fs.x binaries from Quantum ESPRESSO via the quacc library.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ase.dft.kpoints import bandpath
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from quacc import flow, job
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.recipes.espresso._base import run_and_summarize
from quacc.utils.kpts import convert_pmg_kpts
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import EspressoBandsSchema, Filenames, RunSchema, SourceDirectory


@job
def bands_pw_job(
    atoms: Atoms,
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
    make_bandpath: bool = True,
    line_density: float = 20,
    force_gamma: bool = True,
    test_run: bool = False,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic bandstructure calculation with pw.x.

    !!! Note

        First perform a normal SCF calculation [quacc.recipes.espresso.core.static_job][];
        then use this job if you are interested in calculating only the Kohn-Sham states
        for the given set of k-points

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
    make_bandpath
        If True, it returns the primitive cell for your structure and generates
        the high symmetry k-path using Latmer-Munro approach.
        For more information look at
        [pymatgen.symmetry.bandstructure.HighSymmKpath][]
    line_density
        Density of kpoints along the band path if make_bandpath is True
        For more information [quacc.utils.kpts.convert_pmg_kpts][]
    force_gamma
        Forces gamma-centered k-points when using make_bandpath
        For more information [quacc.utils.kpts.convert_pmg_kpts][]
    test_run
        If True, a test run is performed to check that the calculation input_data is correct or
        to generate some files/info if needed.
    additional_fields
        Additional fields to add to the results dictionary.
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
        "input_data": {"control": {"calculation": "bands", "verbosity": "high"}}
    }
    if make_bandpath:
        structure = AseAtomsAdaptor.get_structure(atoms)
        primitive = SpacegroupAnalyzer(structure).get_primitive_standard_structure()
        atoms = primitive.to_ase_atoms()
        calc_defaults["kpts"] = bandpath(
            convert_pmg_kpts(
                {"line_density": line_density}, atoms, force_gamma=force_gamma
            )[0],
            cell=atoms.get_cell(),
        )

    return run_and_summarize(
        atoms,
        template=EspressoTemplate("pw", test_run=test_run, outdir=prev_outdir),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "pw.x bands"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def bands_pp_job(
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
    test_run: bool = False,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to re-order bands and computes bands-related properties with bands.x.
    This allows one to get the bands structure in a more readable way.

    !!! Note

        This requires a previous [quacc.recipes.espresso.bands.bands_pw_job][] calculation.

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
    test_run
        If True, a test run is performed to check that the calculation input_data is correct or
        to generate some files/info if needed.
    additional_fields
        Additional fields to add to the results dictionary.
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
    return run_and_summarize(
        template=EspressoTemplate("bands", test_run=test_run, outdir=prev_outdir),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        additional_fields={"name": "bands.x post-processing"}
        | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def fermi_surface_job(
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
    test_run: bool = False,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to retrieve the fermi surface with fs.x

    !!! Note

        It requires a previous uniform unshifted k-point grid bands calculation.

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
    test_run
        If True, a test run is performed to check that the calculation input_data is correct or
        to generate some files/info if needed.
    additional_fields
        Additional fields to add to the results dictionary.
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
    return run_and_summarize(
        template=EspressoTemplate("fs", test_run=test_run, outdir=prev_outdir),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        additional_fields={"name": "fs.x fermi_surface"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@flow
def bands_flow(
    atoms: Atoms,
    copy_files: (
        SourceDirectory | list[SourceDirectory] | dict[SourceDirectory, Filenames]
    ),
    run_bands_pp: bool = True,
    run_fermi_surface: bool = False,
    make_bandpath: bool = True,
    line_density: float = 20,
    force_gamma: bool = True,
    job_params: dict[str, Any] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> EspressoBandsSchema:
    """
    Function to compute bands structure and fermi surface using pw.x, bands.x and fs.x.

    Consists of the following steps:

    1. A pw.x non-self consistent calculation
        - name: "bands_pw_job"
        - job : [quacc.recipes.espresso.bands.bands_pw_job][]

    2. A bands.x post-processing calculation
        - name: "bands_pp_job"
        - job : [quacc.recipes.espresso.bands.bands_pp_job][]

    3. A fs.x calculation to obtain the fermi surface
        - name: "fermi_surface_job"
        - job : [quacc.recipes.espresso.bands.fermi_surface_job][]

    Parameters
    ----------
    atoms
        The Atoms object.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    run_bands_pp
        If True, a bands.x post-processing calculation will be carried out.
        This allows to re-order bands and computes band-related properties.
    run_fermi_surface
        If True, a fs.x calculation will be carried out.
        This allows to generate the fermi surface of your structure.
        It requires a uniform unshifted k-point grid bands calculation.
    make_bandpath
        If True, it returns the primitive cell for your structure and generates
        the high symmetry k-path using Latmer-Munro approach.
        For more information look at
        [pymatgen.symmetry.bandstructure.HighSymmKpath][]
    line_density
        Density of kpoints along the band path if make_bandpath is True
        For more information [quacc.utils.kpts.convert_pmg_kpts][]
    force_gamma
        Forces gamma-centered k-points when using make_bandpath
        For more information [quacc.utils.kpts.convert_pmg_kpts][]
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    BandsSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    (bands_pw_job_, bands_pp_job_, fermi_surface_job_) = customize_funcs(
        ["bands_pw_job", "bands_pp_job", "fermi_surface_job"],
        [bands_pw_job, bands_pp_job, fermi_surface_job],
        param_swaps=job_params,
        decorators=job_decorators,
    )

    bands_results = bands_pw_job_(
        atoms,
        copy_files,
        make_bandpath=make_bandpath,
        line_density=line_density,
        force_gamma=force_gamma,
    )
    results = {"bands_pw": bands_results}
    bands_results_dir = bands_results["dir_name"]

    if run_bands_pp:
        bands_pp_results = bands_pp_job_(prev_outdir=bands_results_dir)
        results["bands_pp"] = bands_pp_results

    if run_fermi_surface:
        fermi_results = fermi_surface_job_(prev_outdir=bands_results_dir)
        results["fermi_surface"] = fermi_results

    return results
