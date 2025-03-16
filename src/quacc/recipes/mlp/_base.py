"""Base functions for universal machine-learned interatomic potentials."""

from __future__ import annotations

from functools import lru_cache, wraps
from importlib.util import find_spec
from logging import getLogger
from typing import TYPE_CHECKING

from frozendict import frozendict

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase.calculators.calculator import BaseCalculator


LOGGER = getLogger(__name__)


def freezeargs(func):
    """Convert a mutable dictionary into immutable.
    Useful to make sure dictionary args are compatible with cache
    From https://stackoverflow.com/a/53394430
    """

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
    method: Literal["mace-mp-0", "m3gnet", "chgnet", "sevennet", "orb", "fairchem"],
) -> tuple[BaseCalculator, dict[str, Any]]:
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

    Returns
    -------
    Calculator
        The chosen calculator
    """
    import torch

    if not torch.cuda.is_available():
        LOGGER.warning("CUDA is not available to PyTorch. Calculations will be slow.")

    if method.lower() == "m3gnet":
        import matgl
        from matgl import __version__
        from matgl.ext.ase import PESCalculator

        model = matgl.load_model("M3GNet-MP-2021.2.8-DIRECT-PES")
        calc_defaults = {"potential": model, "stress_weight": 1.0 / 160.21766208}
        calc = PESCalculator

    elif method.lower() == "chgnet":
        from chgnet import __version__
        from chgnet.model.dynamics import CHGNetCalculator

        calc_defaults = {}
        calc = CHGNetCalculator

    elif method.lower() == "mace-mp-0":
        from mace import __version__
        from mace.calculators import mace_mp

        calc_defaults = {"default_dtype": "float64"}
        calc = mace_mp

    elif method.lower() == "sevennet":
        from sevenn import __version__
        from sevenn.sevennet_calculator import SevenNetCalculator

        calc_defaults = {}
        calc = SevenNetCalculator

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

        calc = ORBCalculator
        calc_defaults = {"model": pretrained.orb_v2}

    elif method.lower() == "fairchem":
        from fairchem.core import OCPCalculator, __version__

        calc_defaults = {}
        calc = OCPCalculator

    else:
        raise ValueError(f"Unrecognized {method=}.")

    calc.parameters["version"] = __version__

    return calc, calc_defaults
