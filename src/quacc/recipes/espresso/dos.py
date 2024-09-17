"""DOS/ProjWFC recipes for performing dos calculations"""

from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow, job
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.recipes.espresso._base import run_and_summarize
from quacc.recipes.espresso.core import non_scf_job, static_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from ase.atoms import Atoms

    from quacc.types import (
        EspressoDosSchema,
        EspressoProjwfcSchema,
        Filenames,
        RunSchema,
        SourceDirectory,
    )


@job
def dos_job(
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
    Function to carry out a basic dos.x calculation (density of states).
    It is mainly used to extract the charge density and wavefunction from a previous pw.x calculation.
    It generates the total density of states. For more details, please see
    https://www.quantum-espresso.org/Doc/INPUT_DOS.html

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
        If True, the calculation will be run in test mode. This is useful for quickly
        checking if the calculation will run without errors.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        `ase.io.espresso.write_fortran_namelist` for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    return run_and_summarize(
        template=EspressoTemplate("dos", test_run=test_run, outdir=prev_outdir),
        calc_defaults=None,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "dos.x Density-of-States"}
        | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def projwfc_job(
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
    Function to carry out a basic projwfc.x calculation.
    It is mainly used to extract the charge density and wavefunction from a previous pw.x calculation.
    It can generate partial dos, local dos, spilling parameter and more. Fore more details please see
    https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html

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
        If True, the calculation will be run in test mode. This is useful for quickly
        checking if the calculation will run without errors.
    additional_fields
        Additional fields to add to the results dictionary.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        `ase.io.espresso.write_fortran_namelist` for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    return run_and_summarize(
        template=EspressoTemplate("projwfc", test_run=test_run, outdir=prev_outdir),
        calc_defaults=None,
        calc_swaps=calc_kwargs,
        additional_fields={"name": "projwfc.x Projects-wavefunctions"}
        | (additional_fields or {}),
        copy_files=copy_files,
    )


@flow
def dos_flow(
    atoms: Atoms,
    job_decorators: dict[str, Callable | None] | None = None,
    job_params: dict[str, Any] | None = None,
) -> EspressoDosSchema:
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
        Custom parameters to pass to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    DosSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    job_params = job_params or {}
    default_job_params = {
        "static_job": {
            "kspacing": 0.2,
            "input_data": {"system": {"occupations": "tetrahedra"}},
        },
        "non_scf_job": recursive_dict_merge(
            job_params.get("static_job"),
            {
                "kspacing": 0.01,
                "input_data": {
                    "control": {"calculation": "nscf", "verbosity": "high"},
                    "system": {"occupations": "tetrahedra"},
                },
            },
        ),
    }

    static_job_, non_scf_job_, dos_job_ = customize_funcs(
        ["static_job", "non_scf_job", "dos_job"],
        [static_job, non_scf_job, dos_job],
        param_defaults=default_job_params,
        param_swaps=job_params,
        decorators=job_decorators,
    )

    static_results = static_job_(atoms)
    static_results_dir = static_results["dir_name"]
    non_scf_results = non_scf_job_(atoms, prev_outdir=static_results_dir)
    dos_results = dos_job_(prev_outdir=static_results_dir)

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
) -> EspressoProjwfcSchema:
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
        Custom parameters to pass to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    ProjwfcSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    job_params = job_params or {}
    default_job_params = {
        "static_job": {
            "kspacing": 0.2,
            "input_data": {"system": {"occupations": "tetrahedra"}},
        },
        "non_scf_job": recursive_dict_merge(
            job_params.get("static_job"),
            {
                "kspacing": 0.01,
                "input_data": {
                    "control": {"calculation": "nscf", "verbosity": "high"},
                    "system": {"occupations": "tetrahedra"},
                },
            },
        ),
    }
    static_job_, non_scf_job_, projwfc_job_ = customize_funcs(
        ["static_job", "non_scf_job", "projwfc_job"],
        [static_job, non_scf_job, projwfc_job],
        param_defaults=default_job_params,
        param_swaps=job_params,
        decorators=job_decorators,
    )

    static_results = static_job_(atoms)
    static_results_dir = static_results["dir_name"]
    non_scf_results = non_scf_job_(atoms, prev_outdir=static_results_dir)
    projwfc_results = projwfc_job_(prev_outdir=static_results_dir)

    return {
        "static_job": static_results,
        "non_scf_job": non_scf_results,
        "projwfc_job": projwfc_results,
    }
