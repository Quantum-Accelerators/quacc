from __future__ import annotations

import logging
from importlib.util import find_spec
from typing import TYPE_CHECKING

from ase.atoms import Atoms
from monty.dev import requires

from quacc.atoms.core import copy_atoms

LOGGER = logging.getLogger(__name__)

has_geodesic_interpolate = bool(find_spec("geodesic_interpolate"))

if has_geodesic_interpolate:
    from geodesic_interpolate.geodesic import Geodesic
    from geodesic_interpolate.interpolation import redistribute

if TYPE_CHECKING:
    from typing import Literal


@requires(
    has_geodesic_interpolate,
    "geodesic-interpolate must be installed. Refer to the quacc documentation.",
)
def geodesic_interpolate_wrapper(
    reactant: Atoms,
    product: Atoms,
    n_images: int = 10,
    perform_sweep: bool | Literal["auto"] = "auto",
    redistribute_tol: float = 1e-2,
    smoother_tol: float = 2e-3,
    max_iterations: int = 15,
    max_micro_iterations: int = 20,
    morse_scaling: float = 1.7,
    geometry_friction: float = 1e-2,
    distance_cutoff: float = 3.0,
    sweep_cutoff_size: int = 35,
) -> list[Atoms]:
    """
    Interpolates between two geometries and optimizes the path with the geodesic method.

    Parameters
    ----------
    reactant
        The ASE Atoms object representing the initial geometry.
    product
        The ASE Atoms object representing the final geometry.
    n_images
        Number of images for interpolation. Default is 10.
    perform_sweep
        Whether to sweep across the path optimizing one image at a time.
        Default is to perform sweeping updates if there are more than 35 atoms.
    redistribute_tol
        the value passed to the tol keyword argument of
         geodesic_interpolate.interpolation.redistribute. Default is 1e-2.
    smoother_tol
        the value passed to the tol keyword argument of geodesic_smoother.smooth
        or geodesic_smoother.sweep. Default is 2e-3.
    max_iterations
        Maximum number of minimization iterations. Default is 15.
    max_micro_iterations
        Maximum number of micro iterations for the sweeping algorithm. Default is 20.
    morse_scaling
        Exponential parameter for the Morse potential. Default is 1.7.
    geometry_friction
        Size of friction term used to prevent very large changes in geometry. Default is 1e-2.
    distance_cutoff
        Cut-off value for the distance between a pair of atoms to be included in the coordinate system. Default is 3.0.
    sweep_cutoff_size
        Cut off system size that above which sweep function will be called instead of smooth
        in Geodesic.

    Returns
    -------
    list[Atoms]
        A list of ASE Atoms objects representing the smoothed path between the reactant and product geometries.
    """
    reactant = copy_atoms(reactant)
    product = copy_atoms(product)

    # Read the initial geometries.
    chemical_symbols = reactant.get_chemical_symbols()

    # First redistribute number of images. Perform interpolation if too few and subsampling if too many images are given
    raw_interpolated_positions = redistribute(
        chemical_symbols,
        [reactant.positions, product.positions],
        n_images,
        tol=redistribute_tol,
    )

    # Perform smoothing by minimizing distance in Cartesian coordinates with redundant internal metric
    # to find the appropriate geodesic curve on the hyperspace.
    geodesic_smoother = Geodesic(
        chemical_symbols,
        raw_interpolated_positions,
        morse_scaling,
        threshold=distance_cutoff,
        friction=geometry_friction,
    )
    if perform_sweep == "auto":
        perform_sweep = len(chemical_symbols) > sweep_cutoff_size
    if perform_sweep:
        geodesic_smoother.sweep(
            tol=smoother_tol, max_iter=max_iterations, micro_iter=max_micro_iterations
        )
    else:
        geodesic_smoother.smooth(tol=smoother_tol, max_iter=max_iterations)
    return [
        Atoms(symbols=chemical_symbols, positions=geom)
        for geom in geodesic_smoother.path
    ]
