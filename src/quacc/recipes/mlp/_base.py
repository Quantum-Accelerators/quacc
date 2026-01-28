"""Common utility functions for universal machine-learned interatomic potentials."""

from __future__ import annotations

from pathlib import Path

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


@lru_cache
def _get_omat24_references() -> dict[str, float]:
    """
    Fetch formation energy references for OMAT24-trained models from HuggingFace.

    These references come from https://huggingface.co/facebook/UMA/blob/main/references/form_elem_refs.yaml

    Returns
    -------
    dict[str, float]
        Dictionary mapping element symbols to reference energies (eV/atom).
    """
    import yaml
    from huggingface_hub import hf_hub_download

    LOGGER.info("Downloading OMAT24 formation energy references from HuggingFace...")

    # Download the form_elem_refs.yaml file from HuggingFace
    refs_file = hf_hub_download(
        repo_id="facebook/UMA",
        filename="references/form_elem_refs.yaml",
        repo_type="model",
    )

    # Load and extract the omat references
    with Path.open(refs_file) as f:
        refs_data = yaml.safe_load(f)

    omat_refs = refs_data.get("refs", {}).get("omat", {})

    if not omat_refs:
        raise ValueError("Could not find 'refs.omat' in the downloaded reference file.")

    LOGGER.info(f"Loaded OMAT24 references for {len(omat_refs)} elements.")
    return omat_refs


@lru_cache
def _get_mp20_references() -> dict[str, float]:
    """
    Load formation energy references for MP-20 compatible models.

    These references come from matbench-discovery repository:
    https://github.com/janosh/matbench-discovery

    Returns
    -------
    dict[str, float]
        Dictionary mapping element symbols to reference energies (eV/atom).
    """
    import gzip
    import json
    from pathlib import Path

    LOGGER.info("Loading MP-20 formation energy references from local file...")

    # Load from local gzipped JSON file
    refs_file = (
        Path(__file__).parent
        / "references"
        / "2023-02-07-mp-elemental-reference-entries.json.gz"
    )

    if not refs_file.exists():
        raise FileNotFoundError(
            f"MP-20 reference file not found at {refs_file}. "
            "Please ensure the file is in src/quacc/recipes/mlp/references/"
        )

    # Load the gzipped JSON file
    with gzip.open(refs_file, "rt") as f:
        refs_data = json.load(f)

    # Extract element references based on the expected structure
    # The file should contain element references
    if isinstance(refs_data, dict):
        mp20_refs = refs_data
    else:
        raise ValueError(
            f"Unexpected format in MP-20 reference file: {type(refs_data)}"
        )

    LOGGER.info(f"Loaded MP-20 references for {len(mp20_refs)} elements.")
    return mp20_refs


@freezeargs
@lru_cache
def pick_calculator(
    method: Literal[
        "mace-mp", "m3gnet", "chgnet", "tensornet", "sevennet", "orb", "fairchem"
    ],
    use_formation_energy: bool = False,
    references: Literal["MP20", "OMAT24"] | None = None,
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
    use_formation_energy
        If True, wrap the calculator with FormationEnergyCalculator to compute
        formation energies. Requires fairchem-core package to be installed.
        Supported for all calculator types. Default is False.
    references
        Formation energy references to use. Only used if use_formation_energy=True.
        Options:
        - None: Use built-in references from FormationEnergyCalculator (FAIRChem models only)
        - "OMAT24": Use OMAT24 references from https://huggingface.co/facebook/UMA
        - "MP20": Use MP-20 references from matbench-discovery
        Default is None.
    **calc_kwargs
        Custom kwargs for the underlying calculator. Set a value to
        `quacc.Remove` to remove a pre-existing key entirely.

    Returns
    -------
    BaseCalculator
        The instantiated calculator (optionally wrapped with FormationEnergyCalculator)
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

        orb_model = calc_kwargs.get("model", "orb_v3_conservative_inf_omat")
        orbff = getattr(pretrained, orb_model)()
        calc = ORBCalculator(model=orbff, **calc_kwargs)

    elif method.lower() == "fairchem":
        from fairchem.core import FAIRChemCalculator, __version__

        calc = FAIRChemCalculator.from_model_checkpoint(**calc_kwargs)

    else:
        raise ValueError(f"Unrecognized {method=}.")

    calc.parameters["version"] = __version__

    # Wrap with FormationEnergyCalculator if requested
    if use_formation_energy:
        from fairchem.core.calculate.ase_calculator import FormationEnergyCalculator

        # Determine which reference energies to use
        fe_kwargs = {}

        if references == "OMAT24":
            # Use OMAT24 references from HuggingFace
            fe_kwargs["references"] = _get_omat24_references()
        elif references == "MP20":
            # Use MP-20 references from local file
            fe_kwargs["references"] = _get_mp20_references()
        # If references is None, use built-in references from FormationEnergyCalculator
        # (works for FAIRChem models with task_name specified)

        # Wrap with FormationEnergyCalculator using provided kwargs
        calc = FormationEnergyCalculator(calculator=calc, **fe_kwargs)

    return calc
