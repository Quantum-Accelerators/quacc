from __future__ import annotations

from typing import TYPE_CHECKING

import psutil

from quacc import get_settings

if TYPE_CHECKING:
    from typing import Any

_LABEL = "Gaussian"
_LOG_FILE = f"{_LABEL}.log"


def create_gaussian_defaults(
    xc: str = "wb97xd",
    basis: str = "def2tzvp",
    charge: int = 0,
    spin_multiplicity: int = 1,
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
        "command": f"{settings.GAUSSIAN_CMD} < {_LABEL}.com > {_LOG_FILE}",
        "label": _LABEL,
        "mem": "16GB",
        "chk": "Gaussian.chk",
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
