"""
This module, 'bands.py', contains recipes for performing bands and fermi surface calculations using the
bands.x and fs.x binaries from Quantum ESPRESSO via the quacc library.

"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ase.dft.kpoints import bandpath
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from quacc import job
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.recipes.espresso._base import base_fn
from quacc.utils.kpts import convert_pmg_kpts

if TYPE_CHECKING:
    from typing import Any, TypedDict

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema

    class BandsSchema(TypedDict):
        bands: RunSchema
        bands_pp: RunSchema
        fermi_surface: RunSchema


@job
def bands_job(
    atoms: Atoms,
    prev_dir: str | Path,
    run_bands_pp: bool = True,
    run_fermi_surface: bool = False,
    primitive: bool = True,
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    job_params: dict[str, Any] | None = None,
) -> BandsSchema:
    """
    Function to compute bands structure and fermi surface using pw.x, bands.x and fs.x.
    This is all done in a single job as a multi-step process.

    1. A pw.x non-self consistent calculation
        - name: "bands"

    2. A bands.x post-processing calculation
        - name: "bands_pp"

    3. A fs.x calculation to obtain the fermi surface
        - name: "fermi_surface"

    Parameters
    ----------
    atoms
        The Atoms object.
    prev_dir
        Outdir of the previously ran pw.x calculation. This is used to copy
        the entire tree structure of that directory to the working directory
        of this calculation.
    run_bands_pp
        If True, a bands.x post-processing calculation will be carried out.
        This allows to re-order bands and computes band-related properties.
    run_fermi_surface
        If True, a fs.x calculation will be carried out.
        This allows to generate the fermi surface of your structure.
    primitive
        If True, it gets the primitive cell for your structure and generate
        the high symmetry k-path using Latimer-Munro approach.
        For more information look at
        https://pymatgen.org/pymatgen.symmetry.html#pymatgen.symmetry.bandstructure.HighSymmKpath
    parallel_info
        Dictionary containing information about the parallelization of the
        calculation. See the ASE documentation for more information.
    test_run
        If True, a test run is performed to check that the calculation input_data is correct or
        to generate some files/info if needed.
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictinoary where
        the keys are the names of the jobs and the values are dictionaries of parameters.

    Returns
    -------
    BandsSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    results = {}
    bands_kwargs = job_params.get("bands", {})

    bands_result = _bands(
        atoms=atoms,
        prev_dir=prev_dir,
        primitive=primitive,
        parallel_info=parallel_info,
        test_run=test_run,
        **bands_kwargs,
    )

    if run_bands_pp:
        bands_pp_kwargs = job_params.get("bands_pp", {})
        prev_dir = bands_result["dir_name"]
        bands_pp_results = _bands_pp(
            atoms, prev_dir, parallel_info, test_run, **bands_pp_kwargs
        )

    if run_fermi_surface:
        fermi_kwargs = job_params.get("fermi_surface", {})
        prev_dir = bands_result["dir_name"]
        fermi_results = _fermi_surface(
            atoms, prev_dir, parallel_info, test_run, **fermi_kwargs
        )

    results["bands"] = bands_result
    results["bands_pp"] = bands_pp_results if run_bands_pp else None
    results["fermi_surface"] = fermi_results if run_fermi_surface else None

    return results


def _bands(
    atoms: Atoms,
    prev_dir: str | Path,
    primitive: bool = True,
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic bands calculation with pw.x.

    Parameters
    ----------
    atoms
        The Atoms object.
    prev_dir
        Outdir of the previously ran pw.x calculation. This is used to copy
        the entire tree structure of that directory to the working directory
        of this calculation.
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

    calc_defaults = {
        "input_data": {"control": {"calculation": "bands", "verbosity": "high"}}
    }
    if primitive:
        structure = AseAtomsAdaptor.get_structure(atoms)
        primitive = SpacegroupAnalyzer(structure).get_primitive_standard_structure()
        atoms = AseAtomsAdaptor.get_atoms(primitive)
        calc_defaults["kpts"] = bandpath(
            convert_pmg_kpts({"line_density": 20}, atoms)[0], cell=atoms.get_cell()
        )

    return base_fn(
        atoms,
        template=EspressoTemplate("pw", test_run=test_run),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "pw.x bands"},
        copy_files=prev_dir,
    )


def _bands_pp(
    atoms: Atoms,
    prev_dir: str | Path,
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to re-order bands and computes bands-related properties with bands.x.

    Parameters
    ----------
    atoms
        The Atoms object.
    prev_dir
        Outdir of the previously ran pw.x calculation. This is used to copy
        the entire tree structure of that directory to the working directory
        of this calculation.
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

    calc_defaults = {}

    return base_fn(
        atoms,
        template=EspressoTemplate("bands", test_run=test_run),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "bands.x post-processing"},
        copy_files=prev_dir,
    )


def _fermi_surface(
    atoms: Atoms,
    prev_dir: str | Path,
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to retrieve the fermi surface with fs.x

    Parameters
    ----------
    atoms
        The Atoms object.
    prev_dir
        Outdir of the previously ran pw.x calculation. This is used to copy
        the entire tree structure of that directory to the working directory
        of this calculation.
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

    calc_defaults = {}

    return base_fn(
        atoms,
        template=EspressoTemplate("fs", test_run=test_run),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "fs.x fermi_surface"},
        copy_files=prev_dir,
    )
