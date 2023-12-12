"""Base functions for universal machine-learned interatomic potentials."""
from __future__ import annotations

from functools import lru_cache
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase.calculator.calculator import Calculator


@lru_cache
def _pick_calculator(calculator: str, **kwargs) -> Calculator:
    """
    Adapted from `matcalc.util.get_universal_calculator`.
    """

    if calculator.lower().startswith("m3gnet"):
        import matgl
        from matgl.ext.ase import M3GNetCalculator

        model = matgl.load_model("M3GNet-MP-2021.2.8-DIRECT-PES")
        kwargs.setdefault("stress_weight", 1.0 / 160.21766208)
        calc = M3GNetCalculator(potential=model, **kwargs)

    elif calculator.lower() == "chgnet":
        from chgnet.model.dynamics import CHGNetCalculator

        calc = CHGNetCalculator(**kwargs)

    elif calculator.lower() == "mace":
        from mace.calculators import mace_mp

        calc = mace_mp(**kwargs)

    else:
        raise ValueError(f"Unrecognized {calculator=}.")

    return calc
