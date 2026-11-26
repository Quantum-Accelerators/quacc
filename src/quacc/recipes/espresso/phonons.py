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
from quacc.calculators.espresso.phonons import (
    prepare_phonon_path,
    summarize_phonon_dispersion,
    validate_phonon_dispersion_params,
)
from quacc.calculators.espresso.utils import grid_copy_files, grid_prepare_repr
from quacc.recipes.espresso._base import run_and_summarize
from quacc.utils.dicts import recursive_dict_merge
from quacc.wflow_tools.customizers import customize_jobs

if TYPE_CHECKING:
    from collections import UserDict
    from collections.abc import Callable, Sequence
    from typing import Any

    from ase.atoms import Atoms
    from numpy.typing import ArrayLike

    from quacc.types import (
        CopyFiles,
        EspressoPhononDispersionSchema,
        EspressoPhononDosSchema,
        RunSchema,
        SourceDirectory,
    )


_PHONON_ASR_DEFAULTS: dict[str, dict[str, Any]] = {
    "q2r_job": {"input_data": {"input": {"zasr": "crystal"}}},
    "matdyn_job": {"input_data": {"input": {"asr": "crystal"}}},
}


@job
def phonon_job(
    copy_files: CopyFiles | None = None,
    prev_outdir: SourceDirectory | None = None,
    test_run: bool = False,
    use_phcg: bool = False,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
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
        A list of file-transfer specifications. Each specification contains a `source`
        directory and the `filenames` to copy. Glob patterns are supported.
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
    binary = "phcg" if use_phcg else "ph"
    calc_defaults = {
        "input_data": {
            "inputph": {"tr2_ph": 1e-12, "alpha_mix(1)": 0.1, "verbosity": "high"}
        },
        "qpts": (0, 0, 0),
    }

    return run_and_summarize(
        template=EspressoTemplate(binary, test_run=test_run, outdir=prev_outdir),
        calc_defaults=calc_defaults,
        calc_swaps=calc_kwargs,
        additional_fields={"name": f"{binary}.x Phonon"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def q2r_job(
    copy_files: CopyFiles,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
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
        A list of file-transfer specifications. Each specification contains a `source`
        directory and the `filenames` to copy. Glob patterns are supported.
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
        template=EspressoTemplate("q2r"),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        additional_fields={"name": "q2r.x Phonon"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def matdyn_job(
    copy_files: CopyFiles,
    additional_fields: dict[str, Any] | None = None,
    *,
    qpts: ArrayLike | None = None,
    **calc_kwargs: Any,
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
        A list of file-transfer specifications. Each specification contains a `source`
        directory and the `filenames` to copy. Glob patterns are supported.
    additional_fields
        Additional fields to add to the results dictionary.
    qpts
        Ordered q-points for phonon dispersion (``dos=False``). With
        ``input_data={"input": {"q_in_band_form": True}}``, each row contains
        three endpoint coordinates and an integer interpolation count. A count
        of n gives n intervals to the next endpoint; the final count is ignored.
        A zero count skips interpolation to the next endpoint, preserving a
        disconnected path. Otherwise, each row contains three coordinates and
        is evaluated directly without adding points. Set ``q_in_cryst_coord=True``
        in the input namelist for fractional reciprocal coordinates.

        Alternatively, supply QE's native point-count and q-point lines through
        ``additional_cards`` (a string or list of strings), leaving qpts unset.
        These points are independent of the coarse DFPT grid in ph.x. Do not
        set ``nk1/nk2/nk3`` for dispersion; those describe the DOS grid.
    **calc_kwargs
        Additional keyword arguments to pass to the Espresso calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. See the docstring of
        [quacc.calculators.espresso.espresso.Espresso][] for more information.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.Summarize.run][].
        ``results["matdyn_results"]`` contains ``frequencies`` with shape
        (nqpoints, nmodes) in cm^-1, including negative frequencies for imaginary
        modes, and Cartesian ``qpoints`` in units of 2pi/alat. ``distances`` are
        included when ``flfrq + ".gp"`` exists; QE includes path jumps in these
        cumulative distances. ``has_zstar`` is included when the text
        force-constant file is present and records whether it contains the
        dielectric tensor and Born effective charges. The input flags and
        endpoint counts remain in ``parameters``. For DOS, the existing
        ``phonon_dos`` table is retained.

    Examples
    --------
    Evaluate L-Gamma-X-K-Gamma for a primitive fcc cell in fractional reciprocal
    coordinates (81 points):

    ```python
    result = matdyn_job(
        copy_files=[{"source": "q2r_results", "filenames": "q2r.fc*"}],
        input_data={
            "input": {
                "asr": "crystal",
                "dos": False,
                "q_in_band_form": True,
                "q_in_cryst_coord": True,
            }
        },
        qpts=[
            [0.5, 0.5, 0.5, 20],
            [0, 0, 0, 20],
            [0.5, 0, 0.5, 20],
            [0.375, 0.375, 0.75, 20],
            [0, 0, 0, 1],
        ],
    )
    ```

    Set ``q_in_band_form=False`` and pass three-column ``qpts`` to evaluate an
    explicit sequence, including separately sampled disconnected segments.
    For polar insulators, matdyn.x takes the LO-TO direction at Gamma from a
    neighboring q-point; [quacc.recipes.espresso.phonons.phonon_dispersion_flow][]
    arranges Gamma points so that each plotted limit is correct.
    ``flfrq`` and ``fldos`` can be customized in ``input_data["input"]``; their
    defaults remain ``matdyn.freq`` and ``matdyn.dos``. This job sets no acoustic
    sum rule, so matdyn.x applies its own default ``asr="no"`` unless ``asr`` is
    set, as in the example above.
    """
    if qpts is not None:
        calc_kwargs["qpts"] = qpts
    return run_and_summarize(
        template=EspressoTemplate("matdyn"),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        additional_fields={"name": "matdyn.x Phonon"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@flow
def phonon_dos_flow(
    copy_files: CopyFiles | None = None,
    prev_outdir: SourceDirectory | None = None,
    job_params: dict[str, Any] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> EspressoPhononDosSchema:
    """
    Function to carry out a phonon DOS calculation. The phonon calculation is carried
    out on a coarse q-grid, the force constants are calculated and extrapolated to a
    finer q-grid, and the phonon DOS is calculated.

    Between steps, only the files the next program reads are copied: the dynamical
    matrices (`matdyn*`) from ph.x to q2r.x, and the force constants (`q2r.fc*`) from
    q2r.x to matdyn.x.

    By default, the acoustic sum rule `crystal` is applied in q2r.x and matdyn.x. For
    the force constants (matdyn.x `asr`), it imposes the three translational sum rules
    by an optimized projection; for the Born effective charges (q2r.x `zasr`), it
    imposes charge neutrality. This assumes a 3D crystal. Override it through
    `job_params`, for example:

    ```python
    {"matdyn_job": {"input_data": {"input": {"asr": "simple"}}}}
    ```

    - `all` (matdyn.x only): adds rotational invariance and the Huang conditions,
      e.g. for 2D materials or accurate long-wavelength acoustic branches. The
      Huang conditions assume a stress-free structure. For infrared-active solids,
      also set `write_lr=True` in `q2r_job` and `read_lr=True` in `matdyn_job`, e.g.,
      `{"q2r_job": {"input_data": {"input": {"write_lr": True}}}}`.
    - `one-dim` / `zero-dim`: for 1D systems and molecules.

    The Gamma-point frequencies printed by ph.x are not ASR-corrected; a large
    correction may indicate an unconverged calculation.

    Consists of following jobs that can be modified:

    1. ph.x calculation
        - name: "phonon_job"
        - job: [quacc.recipes.espresso.phonons.phonon_job][]
    2. q2r.x calculation
        - name: "q2r_job"
        - job: [quacc.recipes.espresso.phonons.q2r_job][]
        - default: `zasr="crystal"`
    3. matdyn.x calculation
        - name: "matdyn_job"
        - job: [quacc.recipes.espresso.phonons.matdyn_job][]
        - default: `asr="crystal"`

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
        A list of file-transfer specifications. Each specification contains a `source`
        directory and the `filenames` to copy. Glob patterns are supported.
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
    default_job_params: dict[str, dict[str, Any]] = {
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
    default_job_params = recursive_dict_merge(_PHONON_ASR_DEFAULTS, default_job_params)
    ph_job, fc_job, dos_job = customize_jobs(
        {"phonon_job": phonon_job, "q2r_job": q2r_job, "matdyn_job": matdyn_job},
        param_defaults=default_job_params,
        param_swaps=job_params,
        decorators=job_decorators,
    ).values()

    dynamical_matrix_files = "matdyn*"
    real_space_force_constant_files = "q2r.fc*"

    ph_job_results = ph_job(copy_files=copy_files, prev_outdir=prev_outdir)
    fc_job_results = fc_job(
        [{"source": ph_job_results["dir_name"], "filenames": dynamical_matrix_files}]
    )
    dos_job_results = dos_job(
        [
            {
                "source": fc_job_results["dir_name"],
                "filenames": real_space_force_constant_files,
            }
        ]
    )

    return {
        "phonon_job": ph_job_results,
        "q2r_job": fc_job_results,
        "matdyn_job": dos_job_results,
    }


@flow
def phonon_dispersion_flow(
    atoms: Atoms,
    copy_files: CopyFiles | None = None,
    prev_outdir: SourceDirectory | None = None,
    make_bandpath: bool = True,
    line_density: float = 20,
    qpts: ArrayLike | None = None,
    qpoint_labels: Sequence[str] | None = None,
    segment_starts: Sequence[int] | None = None,
    job_params: dict[str, Any] | None = None,
    job_decorators: dict[str, Callable | None] | None = None,
) -> EspressoPhononDispersionSchema:
    """Calculate interpolated harmonic DFPT phonon dispersions from an SCF save.

    Runs the existing ``phonon_job`` (ph.x on a coarse uniform grid),
    ``q2r_job`` (real-space IFCs), and ``matdyn_job`` (path interpolation).
    Only ``matdyn*`` dynamical matrices are transferred to q2r.x and only
    ``q2r.fc*`` force constants to matdyn.x, including compressed/XML files.
    No additional pw.x calculation is run.

    Consists of the following jobs that can be modified:

    1. ph.x calculation
        - name: "phonon_job"
        - job: [quacc.recipes.espresso.phonons.phonon_job][]
        - defaults: `ldisp=True`, `qplot=False`, 4x4x4 q-grid
    2. q2r.x calculation
        - name: "q2r_job"
        - job: [quacc.recipes.espresso.phonons.q2r_job][]
        - default: `zasr="crystal"`
    3. matdyn.x calculation
        - name: "matdyn_job"
        - job: [quacc.recipes.espresso.phonons.matdyn_job][]
        - defaults: `asr="crystal"`, `dos=False`, `q_in_cryst_coord=True`

    Parameters
    ----------
    atoms
        The original SCF structure, with the same cell, orientation, and atom
        count as the converged save. Automatic Latimer-Munro paths use a lattice
        reduction that preserves the cell volume and atom count; coordinates are
        transformed back to this cell's reciprocal basis. No primitive-cell
        conversion is performed. The q-points written by matdyn.x are checked
        against this cell's reciprocal basis.
    copy_files
        File-transfer specifications for the previous SCF calculation, as in
        [quacc.recipes.espresso.phonons.phonon_job][]. Supply this or prev_outdir.
    prev_outdir
        Directory containing the converged SCF save, read directly by ph.x.
        It will be decompressed if necessary.
    make_bandpath
        Generate a high-symmetry path when qpts is unset. If False, qpts is
        required. Supplying qpts takes precedence over automatic generation.
    line_density
        Interpolation intervals per inverse angstrom (including the 2pi factor)
        for automatic paths, following pymatgen's convention. Independent of the
        coarse DFPT grid.
    qpts
        Manual path in fractional reciprocal coordinates of atoms. Each row
        has either three coordinates (an explicit sequence) or coordinates
        plus an integer QE interpolation count (endpoints). A count of n gives
        n intervals to the next endpoint, zero marks a break, and the last
        count is ignored. Endpoint paths are expanded with these counts, and
        matdyn.x always receives an explicit q-point list. For 2D systems in a
        cell with vacuum, an in-plane path can be supplied here.
    qpoint_labels
        One string per manual input point; use empty strings for unlabeled
        points. Automatic paths use pymatgen's Latimer-Munro labels, as in
        [quacc.recipes.espresso.bands.bands_flow][].
    segment_starts
        For explicit three-column qpts, increasing indices beginning with 0
        that mark disconnected segments. For endpoint paths, use zero counts.
        Explicit points are never connected by interpolation.
    job_params
        Overrides for phonon_job, q2r_job, and matdyn_job, including the usual
        "all" key. Set the coarse DFPT mesh through
        ``{"phonon_job": {"input_data": {"inputph": {"nq1": 4, "nq2": 4, "nq3": 4}}}}``.
        The default is 4x4x4, ldisp=True, qplot=False. This mesh must be converged
        for the material; increasing line_density cannot repair unconverged IFCs.
        matdyn uses dos=False and fractional coordinates, without nk1/nk2/nk3.
        Set the path through qpts rather than overriding matdyn's path cards.
    job_decorators
        Custom decorators for the three calculation jobs, including "all".
        Path preparation and result assembly run as standard Quacc jobs, so
        that deferred atoms and results work with workflow engines; they are
        not affected by job_decorators.

    Returns
    -------
    EspressoPhononDispersionSchema
        The three existing job dictionaries plus ``dispersion``: fractional
        qpoints, frequencies of shape (nqpoints, 3*natoms) in cm^-1 (negative
        modes retained), distances in inverse angstroms including the 2pi
        factor, point labels, tick indices/labels/positions, and half-open
        ``segments`` ranges [start, stop). Distances omit inter-segment jumps;
        ticks at a break can share a position while keeping separate labels.
        Metadata records the original cell, atom count, path source, coarse
        grid, ASR, long-range IFC and LO-TO settings (loto_disable, loto_2d),
        and has_zstar, which is true when the force constants contain the
        dielectric tensor and Born effective charges (None if unknown).
        matdyn_job retains QE's raw output for every q-point it evaluated,
        including Cartesian q-points.

    Notes
    -----
    The default q2r zasr and matdyn asr are "crystal", as in phonon_dos_flow,
    and remain configurable. In a dispersion calculation, ph.x computes the
    dielectric tensor and Born effective charges at Gamma automatically when no
    smearing or tetrahedra are used, so matdyn.x applies the non-analytic
    (LO-TO) term. Choose loto_disable or loto_2d explicitly where appropriate;
    loto_2d must be set identically for q2r_job and matdyn_job. The LO-TO term
    is applied at Gamma only when has_zstar is true and neither loto_disable
    nor loto_2d is set; loto_disable changes only the Gamma rows. For long-range
    IFCs, set both q2r write_lr=True and matdyn read_lr=True, as required for
    asr="all" in infrared-active solids. For 2D materials, asr="all" also imposes
    rotational invariance and vanishing stress, which keep the flexural (ZA)
    branch quadratic near Gamma; relax the in-plane lattice first.

    At Gamma, matdyn.x takes the non-analytic direction from the preceding
    q-point, or from the following one when the preceding point is also Gamma.
    Generated paths therefore repeat each interior Gamma point, ending one
    segment with the incoming limit and starting the next with the outgoing
    limit. Where Gamma starts a segment after a break, an auxiliary Gamma point
    is evaluated first and discarded. Manual paths are not split at Gamma; to
    obtain both limits, repeat Gamma and start a new segment at the second copy
    (a zero count on the first copy for endpoint paths, or its index in
    segment_starts for explicit points).
    """
    if copy_files is None and prev_outdir is None:
        raise ValueError(
            "Supply copy_files or prev_outdir from a converged SCF calculation."
        )
    defaults = recursive_dict_merge(
        _PHONON_ASR_DEFAULTS,
        {
            "phonon_job": {
                "input_data": {
                    "inputph": {
                        "ldisp": True,
                        "qplot": False,
                        "nq1": 4,
                        "nq2": 4,
                        "nq3": 4,
                    }
                }
            },
            "matdyn_job": {
                "input_data": {"input": {"dos": False, "q_in_cryst_coord": True}}
            },
        },
    )
    validate_phonon_dispersion_params(defaults, job_params)
    path = job(prepare_phonon_path)(
        atoms, make_bandpath, line_density, qpts, qpoint_labels, segment_starts
    )
    ph_job, fc_job, interpolate_job = customize_jobs(
        {"phonon_job": phonon_job, "q2r_job": q2r_job, "matdyn_job": matdyn_job},
        param_defaults=defaults,
        param_swaps=job_params,
        decorators=job_decorators,
    ).values()
    phonon_results = ph_job(copy_files=copy_files, prev_outdir=prev_outdir)
    q2r_results = fc_job(
        [{"source": phonon_results["dir_name"], "filenames": "matdyn*"}]
    )
    matdyn_results = interpolate_job(
        [{"source": q2r_results["dir_name"], "filenames": "q2r.fc*"}],
        qpts=path["matdyn_qpts"],
    )
    return job(summarize_phonon_dispersion)(
        path, phonon_results, q2r_results, matdyn_results
    )


@flow
def grid_phonon_flow(
    copy_files: CopyFiles | None = None,
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
        A list of file-transfer specifications. Each specification contains a `source`
        directory and the `filenames` to copy. Glob patterns are supported.
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
        prev_dirs = [
            {
                "source": result["dir_name"],
                "filenames": [
                    Path("**", "*.xml.*"),
                    Path("**", "data-file-schema.xml.*"),
                    Path("**", "charge-density.*"),
                    Path("**", "wfc*.*"),
                    Path("**", "paw.txt.*"),
                ],
            }
            for result in grid_results
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
            copy_specs = [
                {"source": source, "filenames": filenames}
                for source, filenames in files_to_copy.items()
            ]
            for representation in repr_to_do:
                ph_input_data["inputph"]["start_irr"] = int(representation[0])
                ph_input_data["inputph"]["last_irr"] = int(representation[-1])
                ph_job_results = ph_job(
                    copy_files=deepcopy(copy_specs), input_data=deepcopy(ph_input_data)
                )
                grid_results.append(ph_job_results)

        return grid_results

    job_params = job_params or {}
    default_job_params: dict[str, dict[str, Any]] = {
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
    ph_init_job, ph_job, ph_recover_job = customize_jobs(
        {"ph_init_job": phonon_job, "ph_job": phonon_job, "ph_recover_job": phonon_job},
        param_defaults=default_job_params,
        param_swaps=job_params,
        decorators=job_decorators,
    ).values()

    ph_init_job_results = ph_init_job(copy_files=copy_files, prev_outdir=prev_outdir)
    grid_results = _grid_phonon_subflow(
        job_params["ph_job"]["input_data"], ph_init_job_results, ph_job, nblocks=nblocks
    )

    return _ph_recover_subflow(grid_results)


@job
def dvscf_q2r_job(
    copy_files: CopyFiles | None = None,
    prev_outdir: SourceDirectory | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
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
        A list of file-transfer specifications. Each specification contains a `source`
        directory and the `filenames` to copy. Glob patterns are supported.
    prev_outdir
        The output directory of a previous calculation. If provided, Quantum Espresso
        will directly read the necessary files from this directory, eliminating the need
        to manually copy files. The directory will be ungzipped if necessary.
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
        template=EspressoTemplate("dvscf_q2r", outdir=prev_outdir),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        additional_fields={"name": "dvscf_q2r Phonon"} | (additional_fields or {}),
        copy_files=copy_files,
    )


@job
def postahc_job(
    copy_files: CopyFiles | None = None,
    prev_outdir: SourceDirectory | None = None,
    additional_fields: dict[str, Any] | None = None,
    **calc_kwargs: Any,
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
        A list of file-transfer specifications. Each specification contains a `source`
        directory and the `filenames` to copy. Glob patterns are supported.
    prev_outdir
        The output directory of a previous calculation. If provided, Quantum Espresso
        will directly read the necessary files from this directory, eliminating the need
        to manually copy files. The directory will be ungzipped if necessary.
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
        template=EspressoTemplate("postahc", outdir=prev_outdir),
        calc_defaults={},
        calc_swaps=calc_kwargs,
        additional_fields={"name": "postahc Phonon"} | (additional_fields or {}),
        copy_files=copy_files,
    )
