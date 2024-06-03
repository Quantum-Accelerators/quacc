"""Transition state recipes for the NewtonNet code."""

from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path
from typing import TYPE_CHECKING

from ase import Atoms
from ase.io import read, write
from ase.mep.neb import NEBOptimizer
from ase.neb import NEB
from monty.dev import requires

from quacc import SETTINGS, change_settings, job, strip_decorator
from quacc.recipes.newtonnet.core import _add_stdev_and_hess, freq_job, relax_job
from quacc.runners.ase import run_opt
from quacc.schemas.ase import summarize_opt_run
from quacc.utils.dicts import recursive_dict_merge

has_geodesic_interpolate = bool(find_spec("geodesic_interpolate"))
has_sella = bool(find_spec("sella"))
has_newtonnet = bool(find_spec("newtonnet"))

if has_sella:
    from sella import IRC, Sella
if has_newtonnet:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet
if has_geodesic_interpolate:
    from geodesic_interpolate.fileio import write_xyz
    from geodesic_interpolate.geodesic import Geodesic
    from geodesic_interpolate.interpolation import redistribute

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.optimize.optimize import Optimizer
    from numpy.typing import NDArray

    from quacc.recipes.newtonnet.core import FreqSchema
    from quacc.runners.ase import OptParams
    from quacc.schemas._aliases.ase import OptSchema

    class TSSchema(OptSchema):
        freq_job: FreqSchema | None

    class IRCSchema(OptSchema):
        freq_job: FreqSchema | None

    class QuasiIRCSchema(OptSchema):
        irc_job: IRCSchema
        freq_job: FreqSchema | None


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
def ts_job(
    atoms: Atoms,
    use_custom_hessian: bool = False,
    run_freq: bool = True,
    freq_job_kwargs: dict[str, Any] | None = None,
    opt_params: OptParams | None = None,
    **calc_kwargs,
) -> TSSchema:
    """
    Perform a transition state (TS) job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    use_custom_hessian
        Whether to use a custom Hessian matrix.
    run_freq
        Whether to run the frequency job.
    freq_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.freq_job][]
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.run_opt][].
    **calc_kwargs
        Dictionary of custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    TSSchema
        Dictionary of results. See the type-hint for the data structure.
    """
    freq_job_kwargs = freq_job_kwargs or {}

    calc_defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
        "hess_method": "autograd",
    }
    opt_defaults = {
        "optimizer": Sella,
        "optimizer_kwargs": (
            {"diag_every_n": 0, "order": 1} if use_custom_hessian else {"order": 1}
        ),
    }

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    atoms.calc = NewtonNet(**calc_flags)

    if use_custom_hessian:
        opt_flags["optimizer_kwargs"]["hessian_function"] = _get_hessian

    atoms.calc = NewtonNet(**calc_flags)

    # Run the TS optimization
    dyn = run_opt(atoms, **opt_flags)
    opt_ts_summary = _add_stdev_and_hess(
        summarize_opt_run(dyn, additional_fields={"name": "NewtonNet TS"})
    )

    # Run a frequency calculation
    freq_summary = (
        strip_decorator(freq_job)(opt_ts_summary["atoms"], **freq_job_kwargs)
        if run_freq
        else None
    )
    opt_ts_summary["freq_job"] = freq_summary

    return opt_ts_summary


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
def irc_job(
    atoms: Atoms,
    direction: Literal["forward", "reverse"] = "forward",
    run_freq: bool = True,
    freq_job_kwargs: dict[str, Any] | None = None,
    opt_params: OptParams | None = None,
    **calc_kwargs,
) -> IRCSchema:
    """
    Perform an intrinsic reaction coordinate (IRC) job using the given atoms object.

    Parameters
    ----------
    atoms
        The atoms object representing the system.
    direction
        The direction of the IRC calculation ("forward" or "reverse").
    run_freq
        Whether to run the frequency analysis.
    freq_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.freq_job][]
    opt_params
        Dictionary of custom kwargs for the optimization process. For a list
        of available keys, refer to [quacc.runners.ase.run_opt][].
    **calc_kwargs
        Custom kwargs for the NewtonNet calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `newtonnet.utils.ase_interface.MLAseCalculator` calculator.

    Returns
    -------
    IRCSchema
        A dictionary containing the IRC summary and thermodynamic summary.
        See the type-hint for the data structure.
    """
    freq_job_kwargs = freq_job_kwargs or {}

    calc_defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    opt_defaults = {
        "optimizer": IRC,
        "optimizer_kwargs": {"dx": 0.1, "eta": 1e-4, "gamma": 0.4, "keep_going": True},
        "run_kwargs": {"direction": direction},
    }

    calc_flags = recursive_dict_merge(calc_defaults, calc_kwargs)
    opt_flags = recursive_dict_merge(opt_defaults, opt_params)

    # Define calculator
    atoms.calc = NewtonNet(**calc_flags)

    # Run IRC
    with change_settings({"CHECK_CONVERGENCE": False}):
        dyn = run_opt(atoms, **opt_flags)
        opt_irc_summary = _add_stdev_and_hess(
            summarize_opt_run(
                dyn, additional_fields={"name": f"NewtonNet IRC: {direction}"}
            )
        )

    # Run frequency job
    freq_summary = (
        strip_decorator(freq_job)(opt_irc_summary["atoms"], **freq_job_kwargs)
        if run_freq
        else None
    )
    opt_irc_summary["freq_job"] = freq_summary

    return opt_irc_summary


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
def quasi_irc_job(
    atoms: Atoms,
    direction: Literal["forward", "reverse"] = "forward",
    run_freq: bool = True,
    irc_job_kwargs: dict[str, Any] | None = None,
    relax_job_kwargs: dict[str, Any] | None = None,
    freq_job_kwargs: dict[str, Any] | None = None,
) -> QuasiIRCSchema:
    """
    Perform a quasi-IRC job using the given atoms object. The initial IRC job by default
    is run with `max_steps: 5`.

    Parameters
    ----------
    atoms
        The atoms object representing the system
    direction
        The direction of the IRC calculation
    run_freq
        Whether to run the frequency analysis
    irc_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.irc_job][]
    relax_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.core.relax_job][]
    freq_job_kwargs
        Keyword arguments to use for the [quacc.recipes.newtonnet.ts.freq_job][]

    Returns
    -------
    QuasiIRCSchema
        A dictionary containing the IRC summary, optimization summary, and
        thermodynamic summary.
        See the type-hint for the data structure.
    """
    relax_job_kwargs = relax_job_kwargs or {}
    freq_job_kwargs = freq_job_kwargs or {}

    irc_job_defaults = {"max_steps": 5}
    irc_job_kwargs = recursive_dict_merge(irc_job_defaults, irc_job_kwargs)

    # Run IRC
    irc_summary = strip_decorator(irc_job)(
        atoms, direction=direction, run_freq=False, **irc_job_kwargs
    )

    # Run opt
    relax_summary = strip_decorator(relax_job)(irc_summary["atoms"], **relax_job_kwargs)

    # Run frequency
    freq_summary = (
        strip_decorator(freq_job)(relax_summary["atoms"], **freq_job_kwargs)
        if run_freq
        else None
    )
    relax_summary["freq_job"] = freq_summary
    relax_summary["irc_job"] = irc_summary

    return relax_summary


def _get_hessian(atoms: Atoms) -> NDArray:
    """
    Calculate and retrieve the Hessian matrix for the given molecular configuration.

    This function takes an ASE Atoms object representing a molecular
    configuration and uses the NewtonNet machine learning calculator to
    calculate the Hessian matrix. The calculated Hessian matrix is then
    returned.

    Parameters
    ----------
    atoms
        The ASE Atoms object representing the molecular configuration.

    Returns
    -------
    NDArray
        The calculated Hessian matrix, reshaped into a 2D array.
    """
    ml_calculator = NewtonNet(
        model_path=SETTINGS.NEWTONNET_MODEL_PATH,
        settings_path=SETTINGS.NEWTONNET_CONFIG_PATH,
    )
    ml_calculator.calculate(atoms)

    return ml_calculator.results["hessian"].reshape((-1, 3 * len(atoms)))


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
@requires(
    has_geodesic_interpolate,
    "geodesic_interpolate must be installed. "
    "git clone https://github.com/virtualzx-nad/geodesic-interpolate.git.",
)
def geodesic_interpolate_wrapper(
    reactant_product_atoms: list[Atoms],
    nimages: int = 20,
    perform_sweep: bool | None = None,
    output_filepath: str | Path = "interpolated.xyz",
    convergence_tolerance: float = 2e-3,
    max_iterations: int = 15,
    max_micro_iterations: int = 20,
    morse_scaling: float = 1.7,
    geometry_friction: float = 1e-2,
    distance_cutoff: float = 3.0,
    save_raw_path: str | Path | None = None,
) -> tuple[list[str], list[list[float]]]:
    """
    Interpolates between two geometries and optimizes the path.

    Parameters:
    -----------
    reactant_product_atoms : List[Atoms]
        List of ASE Atoms objects containing initial and final geometries.
    nimages : int, optional
        Number of images for interpolation. Default is 20.
    perform_sweep : Optional[bool], optional
        Whether to sweep across the path optimizing one image at a time.
        Default is to perform sweeping updates if there are more than 35 atoms.
    output_filepath : Union[str, Path], optional
        Output filename. Default is "interpolated.xyz".
    convergence_tolerance : float, optional
        Convergence tolerance. Default is 2e-3.
    max_iterations : int, optional
        Maximum number of minimization iterations. Default is 15.
    max_micro_iterations : int, optional
        Maximum number of micro iterations for the sweeping algorithm. Default is 20.
    morse_scaling : float, optional
        Exponential parameter for the Morse potential. Default is 1.7.
    geometry_friction : float, optional
        Size of friction term used to prevent very large changes in geometry. Default is 1e-2.
    distance_cutoff : float, optional
        Cut-off value for the distance between a pair of atoms to be included in the coordinate system. Default is 3.0.
    save_raw_path : Optional[Union[str, Path]], optional
        When specified, save the raw path after bisections but before smoothing. Default is None.

    Returns:
    --------
    Tuple[List[str], List[List[float]]]
        A tuple containing the list of symbols and the smoothed path.
    """
    if len(reactant_product_atoms) < 2:
        raise ValueError("Need at least two initial geometries.")

    # Read the initial geometries.
    chemical_symbols = reactant_product_atoms[0].get_chemical_symbols()
    initial_positions = [
        configuration.get_positions() for configuration in reactant_product_atoms
    ]

    # First redistribute number of images. Perform interpolation if too few and subsampling if too many images are given
    raw_interpolated_positions = redistribute(
        chemical_symbols, initial_positions, nimages, tol=convergence_tolerance * 5
    )

    if save_raw_path is not None:
        write_xyz(save_raw_path, chemical_symbols, raw_interpolated_positions)

    # Perform smoothing by minimizing distance in Cartesian coordinates with redundant internal metric
    # to find the appropriate geodesic curve on the hyperspace.
    geodesic_smoother = Geodesic(
        chemical_symbols,
        raw_interpolated_positions,
        morse_scaling,
        threshold=distance_cutoff,
        friction=geometry_friction,
    )
    if perform_sweep is None:
        perform_sweep = len(chemical_symbols) > 35
    try:
        if perform_sweep:
            geodesic_smoother.sweep(
                tol=convergence_tolerance,
                max_iter=max_iterations,
                micro_iter=max_micro_iterations,
            )
        else:
            geodesic_smoother.smooth(tol=convergence_tolerance, max_iter=max_iterations)
    finally:
        # Save the smoothed path to output file. try block is to ensure output is saved if one ^C the process, or there is an error
        write_xyz(output_filepath, chemical_symbols, geodesic_smoother.path)
    return chemical_symbols, geodesic_smoother.path


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
@requires(
    has_geodesic_interpolate,
    "geodesic_interpolate must be installed. "
    "git clone https://github.com/virtualzx-nad/geodesic-interpolate.git.",
)
def setup_images(logdir: str, xyz_r_p: str, n_intermediate: int = 40):
    """
    Sets up intermediate images for NEB calculations between reactant and product states.

    Parameters:
    logdir (str): Directory to save the intermediate files.
    xyz_r_p (str): Path to the XYZ file containing reactant and product structures.
    n_intermediate (int): Number of intermediate images to generate.

    Returns:
    List: List of ASE Atoms objects with calculated energies and forces.
    """
    calc_defaults = {
        "model_path": SETTINGS.NEWTONNET_MODEL_PATH,
        "settings_path": SETTINGS.NEWTONNET_CONFIG_PATH,
    }
    print('calc_defaults:\n\n\n\n\n\n\n\n\n\n\n\n\n', calc_defaults)
    opt_defaults = {"optimizer": Sella, "optimizer_kwargs": ({"order": 0})}
    calc_flags = recursive_dict_merge(calc_defaults, {})
    opt_flags = recursive_dict_merge(opt_defaults, {})

    # try:
    # Ensure the log directory exists
    if logdir is not None:
        Path(logdir).mkdir(parents=True, exist_ok=True)

    # Read reactant and product structures
    reactant = read(xyz_r_p, index="0")
    product = read(xyz_r_p, index="1")

    # Optimize reactant and product structures using sella
    for atom, _name in zip([reactant, product], ["reactant", "product"]):
        # atom.calc = calc()
        atom.calc = NewtonNet(**calc_flags)

        # Run the TS optimization
        dyn = run_opt(atom, **opt_flags)
        opt_ts_summary = _add_stdev_and_hess(
            summarize_opt_run(dyn, additional_fields={"name": "NewtonNet TS"})
        )
        reactant = opt_ts_summary["atoms"].copy()
        # traj_file = Path(logdir) / f"{name}_opt.traj"
        # sella_wrapper(atom, traj_file=traj_file, sella_order=0)
    # Save optimized reactant and product structures
    if logdir is not None:
        r_p_path = Path(logdir) / "r_p.xyz"
        write(r_p_path, [reactant.copy(), product.copy()])

    # Generate intermediate images using geodesic interpolation
    symbols, smoother_path = geodesic_interpolate_wrapper(
        [reactant.copy(), product.copy()], nimages=n_intermediate
    )
    images = [Atoms(symbols=symbols, positions=conf) for conf in smoother_path]

    # Calculate energies and forces for each intermediate image
    for image in images:
        # image.calc = calc()
        # ml_calculator = calc()
        image.calc = NewtonNet(**calc_flags)
        ml_calculator = NewtonNet(**calc_flags)
        ml_calculator.calculate(image)

        energy = ml_calculator.results["energy"]
        forces = ml_calculator.results["forces"]

        image.info["energy"] = energy
        image.arrays["forces"] = forces

    # Save the geodesic path
    if logdir is not None:
        geodesic_path = Path(logdir) / "geodesic_path.xyz"
        write(geodesic_path, images)

    return images

    # except Exception:
    #     return []


@job
@requires(
    has_newtonnet, "NewtonNet must be installed. Refer to the quacc documentation."
)
@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
@requires(
    has_geodesic_interpolate,
    "geodesic_interpolate must be installed. "
    "git clone https://github.com/virtualzx-nad/geodesic-interpolate.git.",
)
def run_neb_method(
    method: str,
    optimizer: type[Optimizer] | None = NEBOptimizer,
    opt_method: str | None = "aseneb",
    precon: str | None = None,
    logdir: str | None = None,
    xyz_r_p: str | None = None,
    n_intermediate: int | None = 20,
    k: float | None = 0.1,
    max_steps: int | None = 1000,
    fmax_cutoff: float | None = 1e-2,
) -> list[Atoms]:
    """
    Run NEB method.

    Args:
        method (str): NEB method.
        optimizer (Optional[Type[Optimizer]]): NEB path optimizer class. Defaults to NEBOptimizer.
        opt_method (Optional[str]): Optimization method. Defaults to 'aseneb'.
        precon (Optional[str]): Preconditioner method. Defaults to None.
        logdir (Optional[str]): Directory to save logs. Defaults to None.
        xyz_r_p (Optional[str]): Path to reactant and product XYZ files. Defaults to None.
        n_intermediate (Optional[int]): Number of intermediate images. Defaults to 20.
        k (Optional[float]): Force constant for the springs in NEB. Defaults to 0.1.
        max_steps (Optional[int]): Maximum number of optimization steps allowed. Defaults to 1000.
        fmax_cutoff (Optional[float]): Convergence cut-off criteria for the NEB optimization. Defaults to 1e-2.

    Returns:
        List[Atoms]: The optimized images.
    """
    images = setup_images(logdir, xyz_r_p, n_intermediate=n_intermediate)

    mep = NEB(
        images,
        k=k,
        method=method,
        climb=True,
        precon=precon,
        remove_rotation_and_translation=True,
        parallel=True,
    )

    if logdir is not None:
        Path(logdir).mkdir(parents=True, exist_ok=True)
        log_filename = f"neb_band_{method}_{optimizer.__name__}_{precon}.txt"
        logfile_path = Path(logdir) / log_filename
    else:
        logfile_path = None

    opt = optimizer(mep, method=opt_method, logfile=logfile_path, verbose=2)

    opt.run(fmax=fmax_cutoff, steps=max_steps)

    # The following was written because of some error in writing the xyz file below
    images_copy = []
    for image in images:
        image_copy = Atoms(symbols=image.symbols, positions=image.positions)
        image_copy.info["energy"] = image.get_potential_energy()
        images_copy.append(image_copy)

    if logdir is not None:
        write(
            f"{logdir}/optimized_path_{method}_{optimizer.__name__}_{precon}.xyz",
            images_copy,
        )

    return images
