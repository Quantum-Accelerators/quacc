"""Base functions for universal machine-learned interatomic potentials."""
from __future__ import annotations

import warnings
from functools import lru_cache
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

    from ase.calculator.calculator import Calculator


@lru_cache
def _pick_calculator(
    calculator: Literal["mace", "m3gnet", "chgnet"], **kwargs
) -> Calculator:
    """
    Adapted from `matcalc.util.get_universal_calculator`.

    Parameters
    ----------
    calculator
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
        warnings.warn("CUDA is not available to PyTorch. Calculations will be slow.")

    if calculator.lower().startswith("m3gnet"):
        import matgl
        from matgl import __version__
        from matgl.ext.ase import M3GNetCalculator

        model = matgl.load_model("M3GNet-MP-2021.2.8-DIRECT-PES")
        kwargs.setdefault("stress_weight", 1.0 / 160.21766208)
        calc = M3GNetCalculator(potential=model, **kwargs)

    elif calculator.lower() == "chgnet":
        from chgnet import __version__
        from chgnet.model.dynamics import CHGNetCalculator

        calc = CHGNetCalculator(**kwargs)

    elif calculator.lower() == "mace":
        from mace import __version__
        from mace.calculators import mace_mp

        calc = mace_mp(**kwargs)

    else:
        raise ValueError(f"Unrecognized {calculator=}.")

    calc.parameters["version"] = __version__

    return calc
