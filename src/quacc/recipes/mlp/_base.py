"""Base functions for universal machine-learned interatomic potentials."""

from __future__ import annotations

import logging
from functools import lru_cache
from typing import TYPE_CHECKING
from warnings import warn

if TYPE_CHECKING:
    from typing import Literal

    from ase.calculators.calculator import Calculator

logger = logging.getLogger(__name__)


@lru_cache
def pick_calculator(
    method: Literal["mace-mp-0", "m3gnet", "chgnet"], **kwargs
) -> Calculator:
    """
    Adapted from `matcalc.util.get_universal_calculator`.

    .. deprecated:: 0.7.6
            method `mace` will be removed in a later version, it is replaced by 'mace-mp-0'
            which more accurately reflects the nature of the model and allows for versioning
            in the future.

    Parameters
    ----------
    method
        Name of the calculator to use
    **kwargs
        Custom kwargs for the underlying calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `mace.calculators.mace_mp`, `chgnet.model.dynamics.CHGNetCalculator`,
        or `matgl.ext.ase.M3GNetCalculator` calculators.

    Returns
    -------
    Calculator
        The chosen calculator
    """
    import torch

    if not torch.cuda.is_available():
        logger.warning("CUDA is not available to PyTorch. Calculations will be slow.")

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

    elif method.lower() in ["mace-mp-0", "mace"]:
        from mace import __version__
        from mace.calculators import mace_mp

        if method.lower() == "mace":
            warn(
                "'mace' is deprecated and support will be removed. Use 'mace-mp-0' instead!",
                DeprecationWarning,
                stacklevel=3,
            )

        if "default_dtype" not in kwargs:
            kwargs["default_dtype"] = "float64"
        calc = mace_mp(**kwargs)

    else:
        raise ValueError(f"Unrecognized {method=}.")

    calc.parameters["version"] = __version__

    return calc
