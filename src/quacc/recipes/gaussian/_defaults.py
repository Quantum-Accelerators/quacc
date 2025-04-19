from __future__ import annotations

from typing import TYPE_CHECKING

import psutil

from quacc import get_settings

if TYPE_CHECKING:
    from typing import Any

_LABEL: str = "Gaussian"


def create_gaussian_defaults(
    *, xc: str, basis: str, charge: int, spin_multiplicity: int
) -> dict[str, Any]:
    """Create the default calculator kwargs for Gaussian.

    Parameters
    ----------
    xc
        Exchange-correlation functional
    basis
        Basis set
    charge
        Charge of the system
    spin_multiplicity
        Multiplicity of the system

    Returns
    -------
    dict[str, Any]
        Dictionary of default calculator kwargs
    """
    settings = get_settings()

    return {
        "command": f"{settings.GAUSSIAN_CMD} < {_LABEL}.com > {_LABEL}.log",
        "label": _LABEL,
        "chk": f"{_LABEL}.chk",
        "mem": "16GB",
        "nprocshared": psutil.cpu_count(logical=False),
        "xc": xc,
        "basis": basis,
        "charge": charge,
        "mult": spin_multiplicity,
        "pop": "CM5",
        "scf": ["maxcycle=250", "xqc"],
        "integral": "ultrafine",
        "nosymmetry": "",
    }
