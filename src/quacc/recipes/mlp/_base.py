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

        To use `orb` method, `pynanoflann` must be installed. To install `pynanoflann`,
        run `pip install "pynanoflann@git+https://github.com/dwastberg/pynanoflann"`.
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

    if not torch.cuda.is_available():
        LOGGER.warning("CUDA is not available to PyTorch. Calculations will be slow.")

    if method.lower() in ("m3gnet", "chgnet", "tensornet"):
        import matgl
        from matgl import __version__
        from matgl.config import BACKEND

        if BACKEND == "DGL":
            from matgl.ext._ase_dgl import PESCalculator
        else:
            from matgl.ext._ase_pyg import PESCalculator

        if method == "m3gnet":
            model = matgl.load_model("M3GNet-MatPES-PBE-v2025.1-PES")
        elif method == "chgnet":
            model = matgl.load_model("CHGNet-MatPES-PBE-2025.2.10-2.7M-PES")
        elif method == "tensornet":
            model = matgl.load_model("TensorNet-MatPES-PBE-v2025.1-PES")

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
        if not find_spec("pynanoflann"):
            raise ImportError(
                """orb-models requires pynanoflann.
                Install pynanoflann with `pip install "pynanoflann@git+https://github.com/dwastberg/pynanoflann"`.
                """
            )
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
