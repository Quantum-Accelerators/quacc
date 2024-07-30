"""
This module, 'phonons.py', contains recipes for performing phonon calculations using the
ph.x binary from Quantum ESPRESSO via the quacc library.

The recipes provided in this module are jobs and flows that can be used to perform
phonon calculations in different fashion.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING

from ase.io.espresso import Namelist

from quacc import Job, flow, job, subflow
from quacc.calculators.espresso.espresso import EspressoTemplate
from quacc.calculators.espresso.utils import grid_copy_files, grid_prepare_repr
from quacc.recipes.espresso._base import run_and_summarize
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs

if TYPE_CHECKING:
    from collections import UserDict
    from typing import Any, Callable

    from quacc.types import (
        EspressoPhononDosSchema,
        Filenames,
        RunSchema,
        SourceDirectory,
    )


@job
def phonon_job(
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
    test_run: bool = False,
    use_phcg: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic ph.x calculation. It should allow you to
    use all the features of the [ph.x binary](https://www.quantum-espresso.org/Doc/INPUT_PH.html)

    `ph.x` calculates the dynamical matrix at a set of q-points within the Density
    Functional Perturbation Theory (DFPT) framework. The dynamical matrix is used to calculate the phonon frequencies and eigenvectors. Various other properties can be
    calculated using other post-processing tools.

    !!! Note

        Phonon calculations rely on a structure that is tightly converged.
        We suggest running a `relax_job` with the following settings:

        ```python
        inputs_data = {
            "control": {"forc_conv_thr": 5.0e-5},
            "electrons": {"conv_thr": 1e-12},
        }
        ```

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
    use_phcg
        If True, the calculation is performed using the `phcg.x` code which uses a faster algorithm.
        It can be used only if you sample the Brillouin Zone at Gamma and you only need the phonon
        modes at Gamma (molecules typically). It cannot be used with spin-polarization, USPP and PAW.
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
            "inputph": {"tr2_ph": 1e-12, "alpha_mix(1)": 0.1, "verbosity": "high"}
        },
        "qpts": (0, 0, 0),
    }

    binary = "phcg" if use_phcg else "ph"

    return run_and_summarize(
        template=EspressoTemplate(binary, test_run=test_run, outdir=prev_outdir),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": f"{binary}.x Phonon"},
        copy_files=copy_files,
    )


@job
def q2r_job(
    copy_files: (
        SourceDirectory | list[SourceDirectory] | dict[SourceDirectory, Filenames]
    ),
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic q2r.x calculation. It should allow you to
    use all the features of the [q2r.x binary](https://www.quantum-espresso.org/Doc/INPUT_Q2R.html#idm51)

    `q2r.x` reads force constant matrices C(q) produced by the `ph.x` code
    for a grid of q-points and calculates the corresponding set
    of interatomic force constants (IFC), C(R)

    Parameters
    ----------
    copy_files
        Source directory or directories to copy files from. If a `SourceDirectory` or a
        list of `SourceDirectory` is provided, this interface will automatically guess
        which files have to be copied over by looking at the binary and `input_data`.
        If a dict is provided, the mode is manual, keys are source directories and values
        are relative path to files or directories to copy. Glob patterns are supported.
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
        template=EspressoTemplate("q2r"),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        additional_fields={"name": "q2r.x Phonon"},
        copy_files=copy_files,
    )


@job
def matdyn_job(
    copy_files: (
        SourceDirectory | list[SourceDirectory] | dict[SourceDirectory, Filenames]
    ),
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic `matdyn.x` calculation. It should allow you to use
    all the features of the [matdyn.x binary](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html#idm138)

    This program calculates the phonon frequencies for a list of generic
    q vectors starting from the interatomic force constants generated
    from the dynamical matrices as written by DFPT phonon code through
    the program `q2r.x`

    Parameters
    ----------
    copy_files
        Source directory or directories to copy files from. If a `SourceDirectory` or a
        list of `SourceDirectory` is provided, this interface will automatically guess
        which files have to be copied over by looking at the binary and `input_data`.
        If a dict is provided, the mode is manual, keys are source directories and values
        are relative path to files or directories to copy. Glob patterns are supported.
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
        template=EspressoTemplate("matdyn"),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        additional_fields={"name": "matdyn Phonon"},
        copy_files=copy_files,
    )


@flow
def phonon_dos_flow(
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
    job_params: dict[str, Any] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> EspressoPhononDosSchema:
    """
    Function to carry out a phonon DOS calculation. The phonon calculation is carried
    out on a coarse q-grid, the force constants are calculated and extrapolated to a
    finer q-grid, and the phonon DOS is calculated.

    Consists of following jobs that can be modified:

    1. ph.x calculation
        - name: "phonon_job"
        - job: [quacc.recipes.espresso.phonons.phonon_job][]
    2. q2r.x calculation
        - name: "q2r_job"
        - job: [quacc.recipes.espresso.phonons.q2r_job][]
    3. matdyn.x calculation
        - name: "matdyn_job"
        - job: [quacc.recipes.espresso.phonons.matdyn_job][]

    !!! Note

        Phonon calculations rely on a structure that is tightly converged.
        We suggest running a `relax_job` with the following settings:

        ```python
        input_data = {
            "control": {"forc_conv_thr": 5.0e-5},
            "electrons": {"conv_thr": 1e-12},
        }
        ```

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
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictionary where the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """
    default_job_params = {
        "phonon_job": {
            "input_data": {
                "inputph": {
                    "tr2_ph": 1e-12,
                    "alpha_mix(1)": 0.1,
                    "verbosity": "high",
                    "ldisp": True,
                    "nq1": 4,
                    "nq2": 4,
                    "nq3": 4,
                }
            }
        },
        "matdyn_job": {
            "input_data": {"input": {"dos": True, "nk1": 32, "nk2": 32, "nk3": 32}}
        },
    }
    ph_job, fc_job, dos_job = customize_funcs(
        ["phonon_job", "q2r_job", "matdyn_job"],
        [phonon_job, q2r_job, matdyn_job],
        param_defaults=default_job_params,
        param_swaps=job_params,
        decorators=job_decorators,
    )

    ph_job_results = ph_job(copy_files=copy_files, prev_outdir=prev_outdir)
    fc_job_results = fc_job(ph_job_results["dir_name"])
    dos_job_results = dos_job(fc_job_results["dir_name"])

    return {
        "phonon_job": ph_job_results,
        "q2r_job": fc_job_results,
        "matdyn_job": dos_job_results,
    }


@flow
def grid_phonon_flow(
    copy_files: (
        SourceDirectory
        | list[SourceDirectory]
        | dict[SourceDirectory, Filenames]
        | None
    ) = None,
    prev_outdir: SourceDirectory | None = None,
    nblocks: int = 1,
    job_params: dict[str, Any] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> RunSchema:
    """
    This function performs grid parallelization of a ph.x calculation. Each
    representation of each q-point is calculated in a separate job, allowing for
    distributed computation across different machines and times.

    The grid parallelization is a technique to make phonon calculation embarrassingly
    parallel. This function should return similar results to
    [quacc.recipes.espresso.phonons.phonon_job][]. If you don't know about
    grid parallelization please consult the Quantum Espresso user manual and
    exemples.

    This approach requires the data of the pw.x calculation to be copied to each job,
    leading to a total data size on the disk of n*m times the size of the pw.x calculation, where:
    - n is the number of q-points
    - m is the number of representations

    In addition to the data produced by each ph.x calculation. This can
    result in large data sizes for systems with many atoms.

    To mitigate this, an optional "nblocks" argument can be provided. This
    groups multiple representations together in a single job, reducing the
    data size by a factor of nblocks, but also reducing the level of parallelization.
    In the case of nblocks = 0, each job will contain all the representations for each q-point.

    Consists of following jobs that can be modified:

    1. ph.x calculation test_run
        - name: "ph_init_job"
        - job: [quacc.recipes.espresso.phonons.phonon_job][]

    2. (n * m) / nblocks ph.x calculations
        - name: "ph_job"
        - job: [quacc.recipes.espresso.phonons.phonon_job][]

    3. ph.x calculation to gather data and diagonalize each dynamical matrix
        - name: "ph_recover_job"
        - job: [quacc.recipes.espresso.phonons.phonon_job][]

    !!! Note

        Phonon calculations rely on a structure that is tightly converged.
        We suggest running a `relax_job` with the following settings:

        ```python
        inputs_data = {
            "control": {"forc_conv_thr": 5.0e-5},
            "electrons": {"conv_thr": 1e-12},
        }
        ```

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
    nblocks
        The number of representations to group together in a single job.
        This will reduce the amount of data produced by a factor of nblocks.
        If nblocks = 0, each job will contain all the representations for a
        single q-point.
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        See the type-hint for the data structure.
    """

    @subflow
    def _ph_recover_subflow(grid_results: list[RunSchema]) -> RunSchema:
        prev_dirs = {}
        for result in grid_results:
            prev_dirs[result["dir_name"]] = [
                Path("**", "*.xml.*"),
                Path("**", "data-file-schema.xml.*"),
                Path("**", "charge-density.*"),
                Path("**", "wfc*.*"),
                Path("**", "paw.txt.*"),
            ]
        return ph_recover_job(copy_files=prev_dirs)

    @subflow
    def _grid_phonon_subflow(
        ph_input_data: UserDict | None,
        ph_init_job_results: RunSchema,
        ph_job: Job,
        nblocks: int = 1,
    ) -> list[RunSchema]:
        """
        This functions is a subflow used in
        [quacc.recipes.espresso.phonons.grid_phonon_flow][].

        Parameters
        ----------
        ph_input_data
            The input data for the phonon calculation.
        ph_init_job_results
            The results of the phonon 'only_init' job.
        ph_job
            The phonon job to be executed.
        nblocks
            The number of blocks for grouping representations.

        Returns
        -------
        list[RunSchema]
            A list of results from each phonon job.
        """
        ph_input_data = Namelist(ph_input_data)
        ph_input_data.to_nested(binary="ph")

        prev_outdir = ph_init_job_results["parameters"]["input_data"]["inputph"][
            "outdir"
        ]

        grid_results = []
        for qnum, qdata in ph_init_job_results["results"].items():
            ph_input_data["inputph"]["start_q"] = qnum
            ph_input_data["inputph"]["last_q"] = qnum
            repr_to_do = grid_prepare_repr(qdata["representations"], nblocks)
            files_to_copy = grid_copy_files(
                ph_input_data, prev_outdir, qnum, qdata["qpoint"]
            )
            for representation in repr_to_do:
                ph_input_data["inputph"]["start_irr"] = representation[0]
                ph_input_data["inputph"]["last_irr"] = representation[-1]
                ph_job_results = ph_job(
                    copy_files=deepcopy(files_to_copy),
                    input_data=deepcopy(ph_input_data),
                )
                grid_results.append(ph_job_results)

        return grid_results

    job_params = job_params or {}
    default_job_params = {
        "ph_init_job": recursive_dict_merge(
            {"input_data": {"inputph": {"lqdir": True, "only_init": True}}},
            job_params.get("ph_job"),
        ),
        "ph_job": {
            "input_data": {
                "inputph": {"lqdir": True, "low_directory_check": True, "recover": True}
            }
        },
        "ph_recover_job": recursive_dict_merge(
            {"input_data": {"inputph": {"recover": True, "lqdir": True}}},
            job_params.get("ph_job"),
        ),
    }
    ph_init_job, ph_job, ph_recover_job = customize_funcs(
        ["ph_init_job", "ph_job", "ph_recover_job"],
        [phonon_job, phonon_job, phonon_job],
        param_defaults=default_job_params,
        param_swaps=job_params,
        decorators=job_decorators,
    )

    ph_init_job_results = ph_init_job(copy_files=copy_files, prev_outdir=prev_outdir)
    grid_results = _grid_phonon_subflow(
        job_params["ph_job"]["input_data"], ph_init_job_results, ph_job, nblocks=nblocks
    )

    return _ph_recover_subflow(grid_results)


@job
def dvscf_q2r_job(
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
    Function to carry out a basic dvscf_q2r calculation allowing phonon potential
    interpolation from coarse to fine q-point grids using Fourier interpolation.
    It should allow you to use all the features of the dvscf_q2r binary which does
    not have an official documentation.

    To use this, run a [quacc.recipes.espresso.phonons.phonon_job][] on a coarse q-point
    grid, dvscf_q2r.x can then be used to inverse Fourier transform the phonon potentials
    to a real-space supercell, you can later run an additional
    [quacc.recipes.espresso.phonons.phonon_job][] with `ldvscf_interpolation = True`
    to Fourier transform the potentials to desired q points.

    Only one card, &input:

    prefix  : Prepended to input/output filenames, default: 'pwscf'
    outdir  : Directory containing input, output, and scratch files.
              In quacc this is always set to the current working directory.
    fildyn  : File where the dynamical matrix is written.
              In quacc this should always be set to 'matdyn'.
    fildvscf : File where the potential variation is written.
               In quacc this should always be set to 'dvscf'.
               (character, Default: 'dvscf')
    wpot_dir : Directory where the w_pot binary files are written.
               In quacc this is always set to outdir / w_pot
    do_long_range : If .true., subtract the long-range part of the potential
                    before interpolation. Requires epsilon and Born effective
                    charge data in _ph0/prefix.phsave/tensor.xml. default: .false.
    do_charge_neutral : If .true., renormalize phonon potential to impose
                    neutrality of Born effective charges. default: .false.
    verbosity : If 'high', write more information to stdout.

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
    return run_and_summarize(
        template=EspressoTemplate("dvscf_q2r", outdir=prev_outdir),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        additional_fields={"name": "dvscf_q2r Phonon"},
        copy_files=copy_files,
    )


@job
def postahc_job(
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
    Function to carry out a basic postahc calculation. It should allow you to
    use all the features of the [postahc.x binary](https://www.quantum-espresso.org/Doc/INPUT_POSTAHC.html#idm11)

    Calculate the phonon-induced electron self-energy in the full matrix form
    at a given temperature. This requires the results of a previous ph.x calculation
    with `electron_phonon='ahc'`

    self energies calculated and printed by `postahc.x`

    - Total self-energy in the on-shell approximation (OSA)
    - Debye-Waller self-energy in the RIA
    - Total Fan self-energy in the OSA
    - Upper Fan self-energy
    - Lower Fan self-energy in the OSA

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
    return run_and_summarize(
        template=EspressoTemplate("postahc", outdir=prev_outdir),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        additional_fields={"name": "postahc Phonon"},
        copy_files=copy_files,
    )
