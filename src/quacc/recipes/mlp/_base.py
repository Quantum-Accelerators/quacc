"""Base functions for universal machine-learned interatomic potentials."""

from __future__ import annotations

from functools import lru_cache
from importlib.util import find_spec
from logging import getLogger
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

    from ase.calculators.calculator import Calculator

LOGGER = getLogger(__name__)


@lru_cache
def pick_calculator(
    method: Literal["mace-mp-0", "m3gnet", "chgnet", "sevennet", "orb", "fairchem"],
    **kwargs,
) -> Calculator:
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
    **kwargs
        Custom kwargs for the underlying calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `mace.calculators.mace_mp`, `chgnet.model.dynamics.CHGNetCalculator`,
        `matgl.ext.ase.M3GNetCalculator`, `sevenn.sevennet_calculator.SevenNetCalculator`, or
        `orb_models.forcefield.calculator.ORBCalculator` calculators.

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

        calc = SevenNetCalculator(**kwargs)

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

        orb_model = kwargs.get("model", "orb_v2")
        orbff = getattr(pretrained, orb_model)()
        calc = ORBCalculator(model=orbff, **kwargs)

    elif method.lower() == "fairchem":
        from fairchem.core import OCPCalculator, __version__

        model_name = kwargs.pop("model_name", None)
        checkpoint_path = kwargs.pop("checkpoint_path", None)

        if model_name is None and checkpoint_path is None:
            raise ValueError(
                "A valid value for model_name, checkpoint_path must be provided for "
                "FAIRChem calculators!"
            )
        if checkpoint_path is None:
            from fairchem.core.models import model_name_to_local_file

            checkpoint_path = model_name_to_local_file(model_name)

        if isinstance(checkpoint_path, str):
            checkpoint_path = Path(checkpoint_path)

        disable_amp = kwargs.pop("disable_amp", False)

        calc = OCPCalculator(checkpoint_path=checkpoint_path, **kwargs)

        if disable_amp:
            calc.trainer.scaler = None

    else:
        raise ValueError(f"Unrecognized {method=}.")

    calc.parameters["version"] = __version__

    return calc
