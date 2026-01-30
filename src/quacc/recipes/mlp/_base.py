"""Common utility functions for universal machine-learned interatomic potentials."""

from __future__ import annotations

from functools import lru_cache, wraps
from importlib.util import find_spec
from logging import getLogger
from typing import TYPE_CHECKING

from ase.units import GPa as _GPa_to_eV_per_A3
from monty.dev import requires

has_frozen = bool(find_spec("frozendict"))
has_fairchem = bool(find_spec("fairchem"))

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    from ase.calculators.calculator import BaseCalculator

    if has_fairchem:
        from fairchem.core.calculate import InferenceBatcher

LOGGER = getLogger(__name__)


@requires(has_frozen, "frozendict must be installed. Run pip install frozendict.")
def freezeargs(func: Callable) -> Callable:
    """
    Convert a mutable dictionary into immutable.
    Useful to make sure dictionary args are compatible with cache
    From https://stackoverflow.com/a/53394430

    Parameters
    ----------
    func
        Function to be wrapped.

    Returns
    -------
    Callable
        Wrapped function with frozen dictionary arguments.
    """
    from frozendict import frozendict

    @wraps(func)
    def wrapped(*args, **kwargs):
        args = (frozendict(arg) if isinstance(arg, dict) else arg for arg in args)
        kwargs = {
            k: frozendict(v) if isinstance(v, dict) else v for k, v in kwargs.items()
        }
        return func(*args, **kwargs)

    return wrapped


def pick_calculator(
    method: Literal[
        "mace-mp", "m3gnet", "chgnet", "tensornet", "sevennet", "orb", "fairchem"
    ],
    **calc_kwargs,
) -> BaseCalculator:
    """
    Adapted from `matcalc.util.get_universal_calculator`.

    !!! Note

        The `orb_models` are licensed under the APACHE license as found at the following
        link: https://github.com/orbital-materials/orb-models

    Parameters
    ----------
    method
        Name of the calculator to use.
    **calc_kwargs
        Custom kwargs for the underlying calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely.

    Returns
    -------
    BaseCalculator
        The instantiated calculator

    Notes
    -----
    When `predict_unit` is provided for fairchem, caching is bypassed to ensure
    thread-safety for concurrent calculations with InferenceBatcher.
    """
    # Handle fairchem with predict_unit separately (uncached for thread safety)
    if method.lower() == "fairchem" and "predict_unit" in calc_kwargs:
        return _create_fairchem_calculator_with_predict_unit(**calc_kwargs)

    # Use cached version for all other cases
    return _pick_calculator_cached(method, **calc_kwargs)


def _create_fairchem_calculator_with_predict_unit(**calc_kwargs) -> BaseCalculator:
    """
    Create a FAIRChemCalculator with a provided predict_unit.

    This is intentionally NOT cached because each concurrent calculation needs
    its own calculator instance to avoid race conditions when using
    InferenceBatcher for batched inference.
    """
    from fairchem.core import FAIRChemCalculator, __version__

    predict_unit = calc_kwargs.pop("predict_unit")
    task_name = calc_kwargs.pop("task_name", None)
    calc = FAIRChemCalculator(
        predict_unit=predict_unit, task_name=task_name, **calc_kwargs
    )
    calc.parameters["version"] = __version__
    return calc


@freezeargs
@lru_cache
def _pick_calculator_cached(
    method: Literal[
        "mace-mp", "m3gnet", "chgnet", "tensornet", "sevennet", "orb", "fairchem"
    ],
    **calc_kwargs,
) -> BaseCalculator:
    """
    Internal cached version of pick_calculator.

    This function is cached to avoid reloading models for repeated calculations.
    """
    import torch

    if not torch.cuda.is_available():
        LOGGER.warning("CUDA is not available to PyTorch. Calculations will be slow.")

    method = method.lower()

    if "m3gnet" in method or "chgnet" in method or "tensornet" in method:
        import matgl
        from matgl import __version__
        from matgl.ext.ase import PESCalculator

        if method == "m3gnet":
            model = matgl.load_model("M3GNet-MatPES-PBE-v2025.1-PES")
        elif method == "chgnet":
            model = matgl.load_model("CHGNet-MatPES-PBE-2025.2.10-2.7M-PES")
        elif method == "tensornet":
            model = matgl.load_model("TensorNet-MatPES-PBE-v2025.1-PES")
        else:
            model = matgl.load_model(method)

        if "stress_weight" not in calc_kwargs:
            calc_kwargs["stress_weight"] = _GPa_to_eV_per_A3

        calc = PESCalculator(potential=model, **calc_kwargs)

    elif method.lower() == "mace-mp":
        from mace import __version__
        from mace.calculators import mace_mp

        if "default_dtype" not in calc_kwargs:
            calc_kwargs["default_dtype"] = "float64"
        calc = mace_mp(**calc_kwargs)

    elif method.lower() == "sevennet":
        from sevenn import __version__
        from sevenn.sevennet_calculator import SevenNetCalculator

        calc = SevenNetCalculator(**calc_kwargs)

    elif method.lower() == "orb":
        from orb_models import __version__
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.calculator import ORBCalculator

        orb_model = calc_kwargs.get("model", "orb_v2")
        orbff = getattr(pretrained, orb_model)()
        calc = ORBCalculator(model=orbff, **calc_kwargs)

    elif method.lower() == "fairchem":
        from fairchem.core import FAIRChemCalculator, __version__

        calc = FAIRChemCalculator.from_model_checkpoint(**calc_kwargs)

    else:
        raise ValueError(f"Unrecognized {method=}.")

    calc.parameters["version"] = __version__

    return calc


# Single batcher cache - reuses same Ray Serve deployment by swapping checkpoints
_current_batcher: InferenceBatcher | None = None
_current_checkpoint_key: tuple | None = None


@requires(has_fairchem, "fairchem must be installed. Run pip install fairchem-core.")
def get_inference_batcher(
    name_or_path: str,
    inference_settings: str = "default",
    device: str | None = None,
    **batcher_kwargs,
):
    """
    Get or create a cached InferenceBatcher for FAIRChem batched inference.

    The batcher is cached and reused across different model checkpoints by swapping
    out the checkpoint without tearing down the Ray Serve deployment. This avoids
    the overhead of repeatedly creating and destroying Ray Serve actors.

    Parameters
    ----------
    name_or_path
        A model name from fairchem.core.pretrained.available_models or a path
        to the checkpoint file.
    inference_settings
        Settings for inference. Can be "default" (general purpose) or "turbo"
        (optimized for speed but requires fixed atomic composition).
    device
        Optional torch device to load the model onto (e.g., 'cuda', 'cpu').
    **batcher_kwargs
        Additional kwargs to pass to InferenceBatcher. Available options:
        - max_batch_size: Maximum number of atoms in a batch (default: 512)
        - batch_wait_timeout_s: Max time to wait for batch (default: 0.1)
        - num_replicas: Number of Ray Serve replicas (default: 1)
        - concurrency_backend_options: Dict with options like {"max_workers": N}

    Returns
    -------
    InferenceBatcher
        A cached InferenceBatcher instance ready for batched inference.

    Examples
    --------
    >>> batcher = get_inference_batcher("uma-s-1")
    >>> # Use batcher.batch_predict_unit in pick_calculator
    >>> calc = pick_calculator(
    ...     "fairchem", predict_unit=batcher.batch_predict_unit, task_name="omat"
    ... )
    """
    global _current_batcher, _current_checkpoint_key

    from fairchem.core.calculate import InferenceBatcher, pretrained_mlip

    # Build checkpoint key (what changes when model changes)
    checkpoint_key = (name_or_path, inference_settings, device)

    # Build full config key (for batcher creation params)
    (
        name_or_path,
        inference_settings,
        device,
        frozenset(batcher_kwargs.items()) if batcher_kwargs else frozenset(),
    )

    # If batcher exists, check if we need to update checkpoint
    if _current_batcher is not None:
        if _current_checkpoint_key == checkpoint_key:
            # Same checkpoint, reuse batcher
            return _current_batcher
        else:
            # Different checkpoint, update it without shutting down
            # Load the new predict unit
            if name_or_path in pretrained_mlip.available_models:
                new_predict_unit = pretrained_mlip.get_predict_unit(
                    name_or_path, inference_settings=inference_settings, device=device
                )
            else:
                new_predict_unit = pretrained_mlip.load_predict_unit(
                    name_or_path, inference_settings=inference_settings, device=device
                )

            # Update the checkpoint on the existing Ray server
            _current_batcher.update_checkpoint(new_predict_unit)
            _current_checkpoint_key = checkpoint_key
            LOGGER.info(
                f"Updated InferenceBatcher checkpoint to {name_or_path}"
            )
            return _current_batcher

    # Load the predict unit
    if name_or_path in pretrained_mlip.available_models:
        predict_unit = pretrained_mlip.get_predict_unit(
            name_or_path, inference_settings=inference_settings, device=device
        )
    else:
        predict_unit = pretrained_mlip.load_predict_unit(
            name_or_path, inference_settings=inference_settings, device=device
        )

    # Get batcher settings with defaults (make a copy to avoid mutating the cache key)
    batcher_kwargs_copy = dict(batcher_kwargs) if batcher_kwargs else {}
    max_batch_size = batcher_kwargs_copy.pop("max_batch_size", 512)
    batch_wait_timeout_s = batcher_kwargs_copy.pop("batch_wait_timeout_s", 0.1)
    num_replicas = batcher_kwargs_copy.pop("num_replicas", 1)
    concurrency_backend_options = batcher_kwargs_copy.pop(
        "concurrency_backend_options", None
    )

    batcher = InferenceBatcher(
        predict_unit=predict_unit,
        max_batch_size=max_batch_size,
        batch_wait_timeout_s=batch_wait_timeout_s,
        num_replicas=num_replicas,
        concurrency_backend_options=concurrency_backend_options,
        **batcher_kwargs_copy,
    )

    _current_batcher = batcher
    _current_checkpoint_key = checkpoint_key

    LOGGER.info(
        f"Created InferenceBatcher for {name_or_path} "
        f"(max_batch_size={max_batch_size}, num_replicas={num_replicas})"
    )

    return batcher


def shutdown_inference_batchers() -> None:
    """
    Shutdown the cached InferenceBatcher instance and clear the cache.

    This should be called when you're done with batched inference to clean up
    Ray Serve resources.
    """
    global _current_batcher, _current_checkpoint_key

    if _current_batcher is not None:
        _current_batcher.shutdown()
        _current_batcher = None
        _current_checkpoint_key = None

    LOGGER.info("InferenceBatcher shut down and cache cleared.")
