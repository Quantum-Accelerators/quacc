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
from quacc.recipes.espresso._base import base_fn
from quacc.recipes.espresso.core import relax_job
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_funcs, strip_decorator

if TYPE_CHECKING:
    from typing import Any, Callable, TypedDict

    from ase.atoms import Atoms

    from quacc.schemas._aliases.ase import RunSchema
    from quacc.utils.files import Filenames, SourceDirectory

    class PhononDosSchema(TypedDict):
        relax_job: RunSchema
        phonon_job: RunSchema
        q2r_job: RunSchema
        matdyn_job: RunSchema


@job
def phonon_job(
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames],
    parallel_info: dict[str] | None = None,
    test_run: bool = False,
    **calc_kwargs,
) -> RunSchema:
    """
    Function to carry out a basic ph.x calculation. It should allow you to
    use all the features of the [ph.x binary](https://www.quantum-espresso.org/Doc/INPUT_PH.html)

    This job requires the results of a previous pw.x calculation, you might
    want to create your own flow to run both jobs in sequence.

    Parameters
    ----------
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
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
        "input_data": {
            "inputph": {"tr2_ph": 1e-12, "alpha_mix(1)": 0.1, "verbosity": "high"}
        },
        "qpts": (0, 0, 0),
    }

    return base_fn(
        template=EspressoTemplate("ph", test_run=test_run),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "ph.x Phonon"},
        copy_files=copy_files,
    )


@job
def q2r_job(
    prev_dir: SourceDirectory, parallel_info: dict[str] | None = None, **calc_kwargs
) -> RunSchema:
    """
    Function to carry out a basic q2r.x calculation. It should allow you to
    use all the features of the [q2r.x binary](https://www.quantum-espresso.org/Doc/INPUT_Q2R.html#idm51)

    This job requires the results of a previous ph.x calculation, you might
    want to create your own flow to run both jobs in sequence.

    Parameters
    ----------
    prev_dir
        Outdir of the previously ran ph.x calculation. This is used to copy
        the the dynamical matrix files.
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

    input_data = Namelist(calc_kwargs.get("input_data"))
    input_data.to_nested(binary="q2r")

    fildyn = input_data["input"].get("fildyn", "matdyn")

    calc_defaults = {"input_data": {"input": {"flfrc": "q2r.fc", "fildyn": fildyn}}}

    copy_files = {prev_dir: [f"{fildyn}*"]}

    return base_fn(
        template=EspressoTemplate("q2r"),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "q2r.x Phonon"},
        copy_files=copy_files,
    )


@job
def matdyn_job(
    prev_dir: SourceDirectory, parallel_info: dict[str] | None = None, **calc_kwargs
) -> RunSchema:
    """
    Function to carry out a basic matdyn.x calculation. It should allow you to use
    all the features of the [matdyn.x binary](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html#idm138)

    This job requires the results of a previous q2r.x calculation, you might
    want to create your own flow to run both jobs in sequence.

    Parameters
    ----------
    prev_dir
        Outdir of the previously ran q2r.x calculation. This is used to copy
        the the force constant file.
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

    input_data = Namelist(calc_kwargs.get("input_data"))
    input_data.to_nested(binary="matdyn")

    flfrc = input_data["input"].get("flfrc", "q2r.fc")

    calc_defaults = {"input_data": {"input": {"flfrc": flfrc}}}

    copy_files = {prev_dir: [f"{flfrc}*"]}

    return base_fn(
        template=EspressoTemplate("matdyn"),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        parallel_info=parallel_info,
        additional_fields={"name": "matdyn Phonon"},
        copy_files=copy_files,
    )


@flow
def phonon_dos_flow(
    atoms: Atoms,
    job_params: dict[str, Any] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> PhononDosSchema:
    """
    Function to carry out a phonon DOS calculation. The phonon calculation is carried out on a coarse q-grid, the force constants are calculated
    and extrapolated to a finer q-grid, and the phonon DOS is calculated.

    Consists of following jobs that can be modified:

    1. pw.x relaxation
        - name: "relax_job"
        - job: [quacc.recipes.espresso.core.relax_job][]
    2. ph.x calculation
        - name: "phonon_job"
        - job: [quacc.recipes.espresso.phonons.phonon_job][]
    3. q2r.x calculation
        - name: "q2r_job"
        - job: [quacc.recipes.espresso.phonons.q2r_job][]
    4. matdyn.x calculation
        - name: "matdyn_job"
        - job: [quacc.recipes.espresso.phonons.matdyn_job][]

    Parameters
    ----------
    atoms
        Atoms object to calculate the phonon DOS.
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictionary where the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    relax_job_defaults = {
        "input_data": {
            "control": {"forc_conv_thr": 5.0e-5},
            "electrons": {"conv_thr": 1e-12},
        }
    }
    ph_job_defaults = {
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
    }
    matdyn_job_defaults = {
        "input_data": {"input": {"dos": True, "nk1": 32, "nk2": 32, "nk3": 32}}
    }

    calc_defaults = {
        "relax_job": relax_job_defaults,
        "phonon_job": ph_job_defaults,
        "matdyn_job": matdyn_job_defaults,
    }

    job_params = recursive_dict_merge(calc_defaults, job_params)

    pw_job, ph_job, fc_job, dos_job = customize_funcs(
        ["relax_job", "phonon_job", "q2r_job", "matdyn_job"],
        [relax_job, phonon_job, q2r_job, matdyn_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    pw_job_results = pw_job(atoms)
    ph_job_results = ph_job(pw_job_results["dir_name"])
    fc_job_results = fc_job(ph_job_results["dir_name"])
    dos_job_results = dos_job(fc_job_results["dir_name"])

    return {
        "relax_job": pw_job_results,
        "phonon_job": ph_job_results,
        "q2r_job": fc_job_results,
        "matdyn_job": dos_job_results,
    }


@flow
def grid_phonon_flow(
    atoms: Atoms,
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

    WARNING: Using the ph.x gamma trick is only partially supported by this function.
    The gamma trick will lead to explicit calculations for each mode. If some of them
    can be calculated using symmetry a full job will still be dispatched. This is
    will become a problem if you system is large: you will dispatch large HPC calculations
    for nothing. This can be tempered by setting a large or zero nblocks value.

    Consists of following jobs that can be modified:

    1. pw.x relaxation
        - name: "relax_job"
        - job: [quacc.recipes.espresso.core.relax_job][]

    2. ph.x calculation test_run
        - name: "ph_init_job"
        - job: [quacc.recipes.espresso.phonons.phonon_job][]

    3. (n * m) / nblocks ph.x calculations
        - name: "ph_job"
        - job: [quacc.recipes.espresso.phonons.phonon_job][]

    4. ph.x calculation to gather data and diagonalize each dynamical matrix
        - name: "ph_recover_job"
        - job: [quacc.recipes.espresso.phonons.phonon_job][]

    Parameters
    ----------
    atoms
        Atoms object
    nblocks
        The number of representations to group together in a single job.
        This will reduce the amount of data produced by a factor of nblocks.
        If nblocks = 0, each job will contain all the representations for a
        single q-point.
    job_params
        Custom parameters to pass to each Job in the Flow. This is a dictinoary where
        the keys are the names of the jobs and the values are dictionaries of parameters.
    job_decorators
        Custom decorators to apply to each Job in the Flow. This is a dictionary where
        the keys are the names of the jobs and the values are decorators.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][].
        See the type-hint for the data structure.
    """

    @job
    def _ph_recover_job(grid_results: list[RunSchema]) -> RunSchema:
        prev_dirs = {}
        for result in grid_results:
            prev_dirs[result["dir_name"]] = [
                Path("**", "*.xml.*"),
                Path("**", "data-file-schema.xml.*"),
                Path("**", "charge-density.*"),
                Path("**", "wfc*.*"),
                Path("**", "paw.txt.*"),
            ]
        return strip_decorator(ph_recover_job)(prev_dirs)

    @subflow
    def _grid_phonon_subflow(
        ph_input_data: dict | None,
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
            The number of blocks for grouping representations. Defaults to 1.

        Returns
        -------
        list[RunSchema]
            A list of results from each phonon job.
        """

        ph_input_data = Namelist(ph_input_data)
        ph_input_data.to_nested(binary="ph")

        grid_results = []
        for qnum, qdata in ph_init_job_results["results"].items():
            ph_input_data["inputph"]["start_q"] = qnum
            ph_input_data["inputph"]["last_q"] = qnum
            repr_to_do = grid_prepare_repr(qdata["representations"], nblocks)
            files_to_copy = grid_copy_files(
                ph_input_data, ph_init_job_results["dir_name"], qnum, qdata["qpoint"]
            )
            for representation in repr_to_do:
                ph_input_data["inputph"]["start_irr"] = representation[0]
                ph_input_data["inputph"]["last_irr"] = representation[-1]
                ph_job_results = ph_job(
                    deepcopy(files_to_copy), input_data=deepcopy(ph_input_data)
                )
                grid_results.append(ph_job_results)

        return grid_results

    job_params = job_params or {}
    relax_job_defaults = {
        "input_data": {
            "control": {"forc_conv_thr": 5.0e-5},
            "electrons": {"conv_thr": 1e-12},
        }
    }
    ph_init_job_defaults = recursive_dict_merge(
        {"input_data": {"inputph": {"lqdir": True, "only_init": True}}},
        job_params.get("ph_job"),
    )
    ph_job_defaults = {
        "input_data": {
            "inputph": {"lqdir": True, "low_directory_check": True, "recover": True}
        }
    }
    ph_recover_job_defaults = recursive_dict_merge(
        {"input_data": {"inputph": {"recover": True, "lqdir": True}}},
        job_params.get("ph_job"),
    )

    calc_defaults = {
        "relax_job": relax_job_defaults,
        "ph_init_job": ph_init_job_defaults,
        "ph_job": ph_job_defaults,
        "ph_recover_job": ph_recover_job_defaults,
    }

    job_params = recursive_dict_merge(calc_defaults, job_params)

    pw_job, ph_init_job, ph_job, ph_recover_job = customize_funcs(
        ["relax_job", "ph_init_job", "ph_job", "ph_recover_job"],
        [relax_job, phonon_job, phonon_job, phonon_job],
        parameters=job_params,
        decorators=job_decorators,
    )

    pw_job_results = pw_job(atoms)

    ph_init_job_results = ph_init_job(pw_job_results["dir_name"])

    grid_results = _grid_phonon_subflow(
        job_params["ph_job"]["input_data"], ph_init_job_results, ph_job, nblocks=nblocks
    )

    return _ph_recover_job(grid_results)
