"""
This module, 'dos.py', contains recipes for performing dos calculations using the
dos.x binary from Quantum ESPRESSO via the quacc library.

The recipes provided in this module are jobs and flows that can be used to perform
dos calculations.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow, job
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.calculators.espresso.utils import pw_copy_files
from quacc.recipes.espresso._base import base_fn
from quacc.recipes.espresso.core import non_scf_job, static_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from typing import Any, Callable, TypedDict

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory

    class DosSchema(TypedDict):
        static_job: RunSchema
        non_scf_job: RunSchema
        dos_job: RunSchema

    class ProjwfcSchema(TypedDict):
        static_job: RunSchema
        non_scf_job: RunSchema


@job
def dos_job(
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames],
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic dos.x calculation (density of states).
    It is mainly used to extract the charge density and wavefunction from a previous pw.x calculation.
    It generates the total density of states. For more details, please see
    https://www.quantum-espresso.org/Doc/INPUT_DOS.html

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
        `ase.io.espresso.write_fortran_namelist` for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    return base_fn(
        template=EspressoTemplate("dos", test_run=test_run),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "dos.x Density-of-States"},
        copy_files=copy_files,
    )


@job
def projwfc_job(
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames],
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic projwfc.x calculation.
    It is mainly used to extract the charge density and wavefunction from a previous pw.x calculation.
    It can generate partial dos, local dos, spilling paramenter and more. Fore more details please see
    https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html

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
        `ase.io.espresso.write_fortran_namelist` for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    return base_fn(
        template=EspressoTemplate("projwfc", test_run=test_run),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "projwfc.x Projects-wavefunctions"},
        copy_files=copy_files,
    )


@flow
def dos_flow(
    atoms: Atoms,
    job_decorators: dict[str, Callable | None] | None = None,
    job_params: dict[str, Any] | None = None,
) -> DosSchema:
    """
    This function performs a total density of states calculations.

    Consists of following jobs that can be modified:

    1. pw.x static
        - name: "static_job"
        - job: [quacc.recipes.espresso.core.static_job][]

    2. pw.x non self-consistent
        - name: "non_scf_job"
        - job: [quacc.recipes.espresso.core.non_scf_job][]

    3. dos.x total density of states
        - name: "dos_job"
        - job: [quacc.recipes.espresso.dos.dos_job][]

    Parameters
    ----------
    atoms
        Atoms object
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictinoary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    DosSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    static_job_defaults = {
        "kspacing": 0.2,
        "input_data": {"system": {"occupations": "tetrahedra"}},
    }
    non_scf_job_defaults = recursive_dict_merge(
        job_params.get("static_job", {}),
        {
            "kspacing": 0.01,
            "input_data": {
                "control": {"calculation": "nscf", "verbosity": "high"},
                "system": {"occupations": "tetrahedra"},
            },
        },
    )
    dos_job_defaults = {}

    calc_defaults = {
        "static_job": static_job_defaults,
        "non_scf_job": non_scf_job_defaults,
        "dos_job": dos_job_defaults,
    }
    job_params = recursive_dict_merge(calc_defaults, job_params)

    static_job_, non_scf_job_, dos_job_ = customize_funcs(
        ["static_job", "non_scf_job", "dos_job"],
        [static_job, non_scf_job, dos_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    static_results = static_job_(atoms)
    files_to_copy = pw_copy_files(
        job_params["static_job"].get("input_data"),
        static_results["dir_name"],
        include_wfc=False,
    )

    non_scf_results = non_scf_job_(atoms, files_to_copy)
    files_to_copy = pw_copy_files(
        job_params["non_scf_job"].get("input_data"),
        non_scf_results["dir_name"],
        include_wfc=False,
    )

    dos_results = dos_job_(files_to_copy)

    return {
        "static_job": static_results,
        "non_scf_job": non_scf_results,
        "dos_job": dos_results,
    }


@flow
def projwfc_flow(
    atoms: Atoms,
    job_decorators: dict[str, Callable | None] | None = None,
    job_params: dict[str, Any] | None = None,
) -> ProjwfcSchema:
    """
    This function performs a projwfc calculation.

    Consists of following jobs that can be modified:

    1. pw.x static
        - name: "static_job"
        - job: [quacc.recipes.espresso.core.static_job][]

    2. pw.x non self-consistent
        - name: "non_scf_job"
        - job: [quacc.recipes.espresso.core.non_scf_job][]

    3. projwfc.x job
        - name: "projwfc_job"
        - job: [quacc.recipes.espresso.dos.projwfc_job][]

    Parameters
    ----------
    atoms
        Atoms object
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictinoary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    ProjwfcSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    static_job_defaults = {
        "kspacing": 0.2,
        "input_data": {"system": {"occupations": "tetrahedra"}},
    }
    non_scf_job_defaults = recursive_dict_merge(
        job_params.get("static_job", {}),
        {
            "kspacing": 0.01,
            "input_data": {
                "control": {"calculation": "nscf", "verbosity": "high"},
                "system": {"occupations": "tetrahedra"},
            },
        },
    )
    projwfc_job_defaults = {}

    calc_defaults = {
        "static_job": static_job_defaults,
        "non_scf_job": non_scf_job_defaults,
        "projwfc_job": projwfc_job_defaults,
    }
    job_params = recursive_dict_merge(calc_defaults, job_params)

    static_job_, non_scf_job_, projwfc_job_ = customize_funcs(
        ["static_job", "non_scf_job", "projwfc_job"],
        [static_job, non_scf_job, projwfc_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    static_results = static_job_(atoms)
    files_to_copy = pw_copy_files(
        job_params["static_job"].get("input_data"),
        static_results["dir_name"],
        include_wfc=False,
    )

    non_scf_results = non_scf_job_(atoms, files_to_copy)
    files_to_copy = pw_copy_files(
        job_params["non_scf_job"].get("input_data"),
        non_scf_results["dir_name"],
        include_wfc=True,
    )

    projwfc_results = projwfc_job_(files_to_copy)

    return {
        "static_job": static_results,
        "non_scf_job": non_scf_results,
        "projwfc_job": projwfc_results,
    }
