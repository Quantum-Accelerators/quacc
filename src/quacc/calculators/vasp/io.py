"""I/O utilities for the Vasp calculator."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from quacc.calculators.vasp import presets
from quacc.utils.files import load_yaml_calc

if TYPE_CHECKING:
    from typing import Any

PRESET_DIR = Path(presets.__file__).parent


def resolve_vasp_yaml_path(name_or_path: str | Path) -> Path:
    """
    Resolve a VASP preset or setups YAML to a file path.

    The lookup order is:

    1. If `name_or_path` is a bare string (no directory component), it is
       matched against the bundled files in `quacc.calculators.vasp.presets`.
       The .yaml/.yml extension is optional, so "DefaultSetGGA" and
       "DefaultSetGGA.yaml" are equivalent.
    2. Otherwise, or if no bundled file matches, it is treated as a path to a
       YAML file on disk. Prefix with "./" or pass a `Path` to bypass the
       bundled lookup.

    Parameters
    ----------
    name_or_path
        The name of a bundled preset or a path to a YAML file.

    Returns
    -------
    Path
        The path to the YAML file.

    Raises
    ------
    FileNotFoundError
        If neither a bundled preset nor an existing file matches.
    """
    if isinstance(name_or_path, str) and Path(name_or_path).name == name_or_path:
        name = Path(name_or_path)
        candidates = (
            [name]
            if name.suffix in (".yaml", ".yml")
            else [Path(f"{name}.yaml"), Path(f"{name}.yml")]
        )
        # Compare against the directory listing so matching is case-sensitive
        # on every OS (Windows/macOS filesystems are case-insensitive)
        bundled = {p.name for p in PRESET_DIR.iterdir() if p.is_file()}
        for candidate in candidates:
            if candidate.name in bundled:
                return PRESET_DIR / candidate

    path = Path(name_or_path).expanduser()
    if path.is_file():
        return path

    available = sorted(p.stem for p in PRESET_DIR.glob("*.y*ml"))
    msg = (
        f"Could not find VASP YAML {str(name_or_path)!r}. It is neither a bundled "
        f"preset nor an existing file. Bundled options: {', '.join(available)}"
    )
    raise FileNotFoundError(msg)


def load_vasp_yaml_calc(yaml_path: str | Path) -> dict[str, Any]:
    """
    Loads a YAML file containing calculator settings. Used for VASP calculations
    and can read quacc-formatted YAMLs that are of the following format:
    ```yaml
    inputs:
      xc: pbe
      algo: all
      setups:
        Cu: Cu_pv
      elemental_magmoms:
        Fe: 5
        Cu: 1
    ```
    where `inputs` is a dictionary of ASE-style input parameters, `setups` is a
    dictionary of ASE-style pseudopotentials, and and `elemental_magmoms` is a
    dictionary of element-wise initial magmoms.

    Parameters
    ----------
    yaml_path
        Path to the YAML file.

    Returns
    -------
    dict
        The calculator configuration (i.e. settings).
    """
    config = load_yaml_calc(yaml_path)

    # Allow for either "Cu_pv" and "_pv" style setups
    if "inputs" in config:
        config["inputs"] = {
            k.lower(): v.lower() if isinstance(v, str) else v
            for k, v in config["inputs"].items()
        }
        setups = config["inputs"].get("setups")
        if isinstance(setups, dict):
            config["inputs"]["setups"] = {
                k: get_setup_suffix(k, v) for k, v in setups.items()
            }

    return config


def get_setup_suffix(symbol: Any, setup: Any) -> Any:
    """
    Convert a full POTCAR name (e.g. "Cu_pv") to the suffix (e.g. "_pv") that
    ASE expects for an element-keyed setup. Suffixes (e.g. "_pv") and setups
    keyed by atom index (which ASE expects as full POTCAR names) are returned
    unchanged.

    Parameters
    ----------
    symbol
        The key of the setups dictionary (an element symbol or atom index).
    setup
        The setup for that key.

    Returns
    -------
    Any
        The setup suffix.
    """
    if isinstance(symbol, str) and isinstance(setup, str):
        return setup.removeprefix(symbol)
    return setup
