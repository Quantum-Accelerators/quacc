"""Utility functions for running ASE calculators with ASE-based methods."""

from __future__ import annotations

import sys
from importlib.util import find_spec
from pathlib import Path
from shutil import copy, copytree
from typing import TYPE_CHECKING, Callable

import numpy as np
from ase import Atoms
from ase.calculators import calculator
from ase.filters import FrechetCellFilter
from ase.io import Trajectory, read, write
from ase.mep import NEB
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from monty.dev import requires
from monty.os.path import zpath

from quacc import SETTINGS
from quacc.atoms.core import copy_atoms, get_final_atoms_from_dynamics
from quacc.runners.prep import calc_cleanup, calc_setup, terminate
from quacc.schemas.ase import summarize_opt_run, summarize_path_opt_run
from quacc.utils.dicts import recursive_dict_merge

has_sella = bool(find_spec("sella"))

if has_sella:
    from sella import Sella

has_newtonnet = bool(find_spec("newtonnet"))

if has_newtonnet:
    from newtonnet.utils.ase_interface import MLAseCalculator as NewtonNet

has_geodesic_interpolate = bool(find_spec("geodesic_interpolate"))

if has_geodesic_interpolate:
    from geodesic_interpolate.fileio import write_xyz
    from geodesic_interpolate.geodesic import Geodesic
    from geodesic_interpolate.interpolation import redistribute

if TYPE_CHECKING:
    from typing import Any, TypedDict

    from ase.optimize.optimize import Optimizer

    from quacc.utils.files import Filenames, SourceDirectory

    class OptParams(TypedDict, total=False):
        """
        Type hint for `opt_params` used throughout quacc.
        """

        fmax: float
        max_steps: int
        optimizer: Optimizer = BFGS
        optimizer_kwargs: OptimizerKwargs | None
        store_intermediate_results: bool
        fn_hook: Callable | None
        run_kwargs: dict[str, Any] | None

    class OptimizerKwargs(TypedDict, total=False):
        """
        Type hint for `optimizer_kwargs` in [quacc.runners.ase.run_opt][].
        """

        restart: Path | str | None  # default = None
        append_trajectory: bool  # default = False

    class VibKwargs(TypedDict, total=False):
        """
        Type hint for `vib_kwargs` in [quacc.runners.ase.run_vib][].
        """

        indices: list[int] | None  # default = None
        delta: float  # default = 0.01
        nfree: int  # default = 2


def run_calc(
    atoms: Atoms,
    geom_file: str | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
    properties: list[str] | None = None,
) -> Atoms:
    """
    Run a calculation in a scratch directory and copy the results back to the original
    directory. This can be useful if file I/O is slow in the working directory, so long
    as file transfer speeds are reasonable.

    This is a wrapper around atoms.get_potential_energy(). Note: This function
    does not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    geom_file
        The filename of the log file that contains the output geometry, used to
        update the atoms object's positions and cell after a job. It is better
        to specify this rather than relying on ASE's
        atoms.get_potential_energy() function to update the positions, as this
        varies between codes.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.
    properties
        List of properties to calculate. Defaults to ["energy"] if `None`.

    Returns
    -------
    Atoms
        The updated Atoms object.
    """
    # Copy atoms so we don't modify it in-place
    atoms = copy_atoms(atoms)

    # Perform staging operations
    tmpdir, job_results_dir = calc_setup(atoms, copy_files=copy_files)

    # Run calculation
    if properties is None:
        properties = ["energy"]

    try:
        atoms.calc.calculate(atoms, properties, calculator.all_changes)
    except Exception as exception:
        terminate(tmpdir, exception)

    # Most ASE calculators do not update the atoms object in-place with a call
    # to .get_potential_energy(), which is important if an internal optimizer is
    # used. This section is done to ensure that the atoms object is updated to
    # the final geometry if `geom_file` is provided.
    # Note: We have to be careful to make sure we don't lose the calculator
    # object, as this contains important information such as the parameters
    # and output properties (e.g. final magnetic moments).
    if geom_file:
        atoms_new = read(zpath(tmpdir / geom_file))
        if isinstance(atoms_new, list):
            atoms_new = atoms_new[-1]

        # Make sure the atom indices didn't get updated somehow (sanity check).
        # If this happens, there is a serious problem.
        if (
            np.array_equal(atoms_new.get_atomic_numbers(), atoms.get_atomic_numbers())
            is False
        ):
            raise ValueError("Atomic numbers do not match between atoms and geom_file.")

        atoms.positions = atoms_new.positions
        atoms.cell = atoms_new.cell

    # Perform cleanup operations
    calc_cleanup(atoms, tmpdir, job_results_dir)

    return atoms


def run_opt(
    atoms: Atoms,
    relax_cell: bool = False,
    fmax: float = 0.01,
    max_steps: int = 1000,
    optimizer: Optimizer = BFGS,
    optimizer_kwargs: OptimizerKwargs | None = None,
    store_intermediate_results: bool = False,
    fn_hook: Callable | None = None,
    run_kwargs: dict[str, Any] | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> Optimizer:
    """
    Run an ASE-based optimization in a scratch directory and copy the results back to
    the original directory. This can be useful if file I/O is slow in the working
    directory, so long as file transfer speeds are reasonable.

    This is a wrapper around the optimizers in ASE. Note: This function does not
    modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    relax_cell
        Whether to relax the unit cell shape and volume.
    fmax
        Tolerance for the force convergence (in eV/A).
    max_steps
        Maximum number of steps to take.
    optimizer
        Optimizer class to use.
    optimizer_kwargs
        Dictionary of kwargs for the optimizer. Takes all valid kwargs for ASE
        Optimizer classes. Refer to `_set_sella_kwargs` for Sella-related
        kwargs and how they are set.
    store_intermediate_results
        Whether to store the files generated at each intermediate step in the
        optimization. If enabled, they will be stored in a directory named
        `stepN` where `N` is the step number, starting at 0.
    fn_hook
        A custom function to call after each step of the optimization.
        The function must take the instantiated dynamics class as
        its only argument.
    run_kwargs
        Dictionary of kwargs for the run() method of the optimizer.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    Optimizer
        The ASE Optimizer object.
    """
    # Copy atoms so we don't modify it in-place
    atoms = copy_atoms(atoms)

    # Perform staging operations
    tmpdir, job_results_dir = calc_setup(atoms, copy_files=copy_files)

    # Set defaults
    optimizer_kwargs = recursive_dict_merge(
        {
            "logfile": "-" if SETTINGS.DEBUG else tmpdir / "opt.log",
            "restart": tmpdir / "opt.json",
        },
        optimizer_kwargs,
    )
    run_kwargs = run_kwargs or {}

    # Check if trajectory kwarg is specified
    if "trajectory" in optimizer_kwargs:
        msg = "Quacc does not support setting the `trajectory` kwarg."
        raise ValueError(msg)

    # Handle optimizer kwargs
    if optimizer.__name__.startswith("SciPy"):
        optimizer_kwargs.pop("restart", None)
    elif optimizer.__name__ == "Sella":
        _set_sella_kwargs(atoms, optimizer_kwargs)
    elif optimizer.__name__ == "IRC":
        optimizer_kwargs.pop("restart", None)

    # Define the Trajectory object
    traj_file = tmpdir / "opt.traj"
    traj = Trajectory(traj_file, "w", atoms=atoms)
    optimizer_kwargs["trajectory"] = traj

    # Set volume relaxation constraints, if relevant
    if relax_cell and atoms.pbc.any():
        atoms = FrechetCellFilter(atoms)

    # Run optimization
    try:
        with traj, optimizer(atoms, **optimizer_kwargs) as dyn:
            if optimizer.__name__.startswith("SciPy"):
                # https://gitlab.com/ase/ase/-/issues/1475
                dyn.run(fmax=fmax, steps=max_steps, **run_kwargs)
            else:
                for i, _ in enumerate(
                    dyn.irun(fmax=fmax, steps=max_steps, **run_kwargs)
                ):
                    if store_intermediate_results:
                        _copy_intermediate_files(
                            tmpdir,
                            i,
                            files_to_ignore=[
                                traj_file,
                                optimizer_kwargs["restart"],
                                optimizer_kwargs["logfile"],
                            ],
                        )
                    if fn_hook:
                        fn_hook(dyn)
    except Exception as exception:
        terminate(tmpdir, exception)

    # Store the trajectory atoms
    dyn.traj_atoms = read(traj_file, index=":")

    # Perform cleanup operations
    calc_cleanup(get_final_atoms_from_dynamics(dyn), tmpdir, job_results_dir)

    return dyn


def run_vib(
    atoms: Atoms,
    vib_kwargs: VibKwargs | None = None,
    copy_files: SourceDirectory | dict[SourceDirectory, Filenames] | None = None,
) -> Vibrations:
    """
    Run an ASE-based vibration analysis in a scratch directory and copy the results back
    to the original directory. This can be useful if file I/O is slow in the working
    directory, so long as file transfer speeds are reasonable.

    This is a wrapper around the vibrations module in ASE. Note: This function
    does not modify the atoms object in-place.

    Parameters
    ----------
    atoms
        The Atoms object to run the calculation on.
    vib_kwargs
        Dictionary of kwargs for the [ase.vibrations.Vibrations][] class.
    copy_files
        Files to copy (and decompress) from source to the runtime directory.

    Returns
    -------
    Vibrations
        The updated Vibrations module
    """
    # Copy atoms so we don't modify it in-place
    atoms = copy_atoms(atoms)

    # Set defaults
    vib_kwargs = vib_kwargs or {}

    # Perform staging operations
    tmpdir, job_results_dir = calc_setup(atoms, copy_files=copy_files)

    # Run calculation
    vib = Vibrations(atoms, name=str(tmpdir / "vib"), **vib_kwargs)
    try:
        vib.run()
    except Exception as exception:
        terminate(tmpdir, exception)

    # Summarize run
    vib.summary(log=sys.stdout if SETTINGS.DEBUG else str(tmpdir / "vib_summary.log"))

    # Perform cleanup operations
    calc_cleanup(vib.atoms, tmpdir, job_results_dir)

    return vib


def run_path_opt(
    xyz_r_p,
    logdir=None,
    optimizer_class=None,
    n_intermediate: int | None = 20,
) -> list[Atoms]:
    """
    Run NEB-based path optimization in a scratch directory and copy the results back to
    the original directory.

    Parameters
    ----------
    reactant_product_path : str
        Path to the XYZ file containing reactant and product structures.
    logdir : str
        Directory to save logs and intermediate files.
    num_intermediate_images : int, optional
        Number of intermediate images to generate. Default is 20.
    spring_constant : float, optional
        Force constant for the springs in NEB. Default is 0.1.
    max_optimization_steps : int, optional
        Maximum number of optimization steps allowed. Default is 1000.
    force_convergence_tolerance : float, optional
        Convergence tolerance for the forces (in eV/A). Default is 0.01.
    neb_method : str, optional
        NEB method to use. Default is 'aseneb'.
    optimizer_class : type[Optimizer], optional
        NEB path optimizer class. Default is NEBOptimizer.
    preconditioner : str | None, optional
        Preconditioner method. Default is None.
    store_intermediate_results : bool, optional
        Whether to store intermediate results at each step. Default is False.
    fn_hook : Callable | None, optional
        Custom function to call after each optimization step. Default is None.
    optimizer_kwargs : dict[str, Any] | None, optional
        Dictionary of kwargs for the optimizer. Default is None.
    run_kwargs : dict[str, Any] | None, optional
        Dictionary of kwargs for the run() method of the optimizer. Default is None.
    copy_files : SourceDirectory | dict[SourceDirectory, Filenames] | None, optional
        Files to copy (and decompress) from source to the runtime directory. Default is None.

    Returns
    -------
    list[Atoms]
        The optimized images.
    """
    # Generate intermediate images
    images = _setup_images(logdir, xyz_r_p, n_intermediate)

    neb = NEB(images)
    neb.interpolate()

    # qn = BFGS(neb, trajectory='neb.traj')
    qn = optimizer_class(neb, trajectory="neb.traj")
    qn.run(fmax=0.05)
    traj = read("neb.traj", ":")

    neb_summary = summarize_path_opt_run(traj, neb, qn)

    return images, neb_summary


def _geodesic_interpolate_wrapper(
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


def _setup_images(logdir: str, xyz_r_p: str, n_intermediate: int = 40):
    """
    Sets up intermediate images for NEB calculations between reactant and product states.

    Parameters:
    logdir (str): Directory to save the intermediate files.
    xyz_r_p (str): Path to the XYZ file containing reactant and product structures.
    n_intermediate (int): Number of intermediate images to generate.

    Returns:
    List: List of ASE Atoms objects with calculated energies and forces.
    """
    current_file_path = Path(__file__).parent
    conf_path = (
        current_file_path / "../../../tests/core/recipes/newtonnet_recipes"
    ).resolve()
    NEWTONNET_CONFIG_PATH = conf_path / "config0.yml"
    NEWTONNET_MODEL_PATH = conf_path / "best_model_state.tar"
    SETTINGS.CHECK_CONVERGENCE = False
    calc_defaults = {
        "model_path": NEWTONNET_MODEL_PATH,
        "settings_path": NEWTONNET_CONFIG_PATH,
    }
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
        atom.calc = NewtonNet(**calc_flags)

        # Run the TS optimization
        dyn = run_opt(atom, **opt_flags)
        opt_ts_summary = summarize_opt_run(
            dyn, additional_fields={"name": "NewtonNet TS"}
        )

        reactant = opt_ts_summary["atoms"].copy()
        # traj_file = Path(logdir) / f"{name}_opt.traj"
        # sella_wrapper(atom, traj_file=traj_file, sella_order=0)
    # Save optimized reactant and product structures
    if logdir is not None:
        r_p_path = Path(logdir) / "r_p.xyz"
        write(r_p_path, [reactant.copy(), product.copy()])

    # Generate intermediate images using geodesic interpolation
    symbols, smoother_path = _geodesic_interpolate_wrapper(
        [reactant.copy(), product.copy()], nimages=n_intermediate
    )
    images = [Atoms(symbols=symbols, positions=conf) for conf in smoother_path]

    # Calculate energies and forces for each intermediate image
    for image in images:
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


@requires(has_sella, "Sella must be installed. Refer to the quacc documentation.")
def _set_sella_kwargs(atoms: Atoms, optimizer_kwargs: dict[str, Any]) -> None:
    """
    Modifies the `optimizer_kwargs` in-place to address various Sella-related
    parameters. This function does the following for the specified key/value pairs in
    `optimizer_kwargs`:

    1. Sets `order = 0` if not specified (i.e. minimization rather than TS
    by default).

    2. If `internal` is not defined and not `atoms.pbc.any()`, set it to `True`.

    Parameters
    ----------
    atoms
        The Atoms object.
    optimizer_kwargs
        The kwargs for the Sella optimizer.

    Returns
    -------
    None
    """
    if "order" not in optimizer_kwargs:
        optimizer_kwargs["order"] = 0

    if not atoms.pbc.any() and "internal" not in optimizer_kwargs:
        optimizer_kwargs["internal"] = True


def _copy_intermediate_files(
    tmpdir: Path, step_number: int, files_to_ignore: list[Path] | None = None
) -> None:
    """
    Copy all files in the working directory to a subdirectory named `stepN` where `N`
    is the step number. This is useful for storing intermediate files generated during
    an ASE relaaxation.

    Parameters
    ----------
    tmpdir
        The working directory.
    step_number
        The step number.
    files_to_ignore
        A list of files to ignore when copying files to the subdirectory.

    Returns
    -------
    None
    """
    files_to_ignore = files_to_ignore or []
    store_path = tmpdir / f"step{step_number}"
    store_path.mkdir()
    for item in tmpdir.iterdir():
        if not item.name.startswith("step") and item not in files_to_ignore:
            if item.is_file():
                copy(item, store_path)
            elif item.is_dir():
                copytree(item, store_path / item.name)
