"""Common utility functions for universal machine-learned interatomic potentials."""

from __future__ import annotations

from functools import lru_cache, wraps
from importlib.util import find_spec
from logging import getLogger
from typing import TYPE_CHECKING

from ase.units import GPa as _GPa_to_eV_per_A3
from monty.dev import requires

has_frozen = bool(find_spec("frozendict"))

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    from ase.calculators.calculator import BaseCalculator

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


@freezeargs
@lru_cache
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
    """
    import torch

    from quacc import get_settings

    settings = get_settings()
    method = method.lower()

    # Skip CUDA warning for rayserve batching mode (inference happens on remote GPU)
    use_ray_serve = method == "fairchem" and settings.FAIRCHEM_RAY_SERVE_BATCHING
    if not use_ray_serve and not torch.cuda.is_available():
        LOGGER.warning("CUDA is not available to PyTorch. Calculations will be slow.")

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

        # Check if Ray Serve batching is enabled AND Ray is actually available
        if use_ray_serve:
            try:
                import ray

                if not ray.is_initialized():
                    LOGGER.warning(
                        "FAIRCHEM_RAY_SERVE_BATCHING is enabled but Ray is not initialized. "
                        "Falling back to local inference. To use Ray Serve batching, ensure "
                        "your flow is running with a Ray cluster (e.g., SlurmRayTaskRunner with "
                        "start_inference_server=True)."
                    )
                    use_ray_serve = False
            except ImportError:
                LOGGER.warning(
                    "FAIRCHEM_RAY_SERVE_BATCHING is enabled but Ray is not installed. "
                    "Falling back to local inference."
                )
                use_ray_serve = False

        if use_ray_serve:
            # Use Ray Serve multiplexed deployment for inference. The deployment
            # is expected to already be running on the cluster (typically started
            # by get_local_ray_cluster / get_slurm_ray_cluster with
            # setup_multiplexed_batch_predict_server). We connect by deployment
            # name and route requests to the appropriate model via the
            # multiplexed_model_id, which has the form
            # "<checkpoint_name_or_path>:<inference_settings>".
            from fairchem.core.units.mlip_unit.predict import BatchServerPredictUnit

            calc_kwargs = dict(calc_kwargs)

            # Determine model identifier: prefer name_or_path (local checkpoint)
            # over model_id/checkpoint
            name_or_path = calc_kwargs.pop("name_or_path", None)
            if name_or_path is not None:
                checkpoint_id = str(name_or_path)
            else:
                checkpoint_id = calc_kwargs.pop("model_id", None) or calc_kwargs.pop(
                    "checkpoint", "uma-s-1p1"
                )

            inference_settings = calc_kwargs.pop("inference_settings", "default")
            task_name = calc_kwargs.pop("task_name", "omat")

            # Drop kwargs only meaningful when loading the checkpoint locally
            calc_kwargs.pop("device", None)
            calc_kwargs.pop("overrides", None)
            calc_kwargs.pop("seed", None)

            multiplexed_model_id = f"{checkpoint_id}:{inference_settings}"

            mlip_unit = BatchServerPredictUnit.from_deployment_connection_info(
                deployment_name="predict-server",
                multiplexed_model_id=multiplexed_model_id,
            )

            calc = FAIRChemCalculator(predict_unit=mlip_unit, task_name=task_name)
        else:
            # Use local inference
            calc = FAIRChemCalculator.from_model_checkpoint(**calc_kwargs)

    else:
        raise ValueError(f"Unrecognized {method=}.")

    calc.parameters["version"] = __version__

    return calc
