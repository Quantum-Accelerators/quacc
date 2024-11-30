"""Base functions for universal machine-learned interatomic potentials."""

from __future__ import annotations

from functools import lru_cache
from logging import getLogger
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

    from ase.calculators.calculator import Calculator

LOGGER = getLogger(__name__)


@lru_cache
def pick_calculator(
    method: Literal["mace-mp-0", "m3gnet", "chgnet", "sevennet"], **kwargs
) -> Calculator:
    """
    Adapted from `matcalc.util.get_universal_calculator`.

    Parameters
    ----------
    method
        Name of the calculator to use
    **kwargs
        Custom kwargs for the underlying calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `mace.calculators.mace_mp`, `chgnet.model.dynamics.CHGNetCalculator`,
        `matgl.ext.ase.M3GNetCalculator`, or `sevenn.sevennet_calculator.SevenNetCalculator`
        calculators.

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
        kwargs.setdefault("stress_weight", 1.0 / 160.21766208)
        calc = PESCalculator(potential=model, **kwargs)

    elif method.lower() == "chgnet":
        from chgnet import __version__
        from chgnet.model.dynamics import CHGNetCalculator

        calc = CHGNetCalculator(**kwargs)

    elif method.lower() == "mace-mp-0":
        from mace import __version__
        from mace.calculators import mace_mp

        if "default_dtype" not in kwargs:
            kwargs["default_dtype"] = "float64"
        calc = mace_mp(**kwargs)
    
    elif method.lower() == "sevennet":
        from sevenn import __version__
        from sevenn.sevennet_calculator import SevenNetCalculator
        
        if "model" not in kwargs:
            kwargs["model"] = "7net-0"
        if "device" not in kwargs:
            kwargs["device"] = "cpu"
        calc = SevenNetCalculator(**kwargs)

    else:
        raise ValueError(f"Unrecognized {method=}.")

    calc.parameters["version"] = __version__

    return calc
