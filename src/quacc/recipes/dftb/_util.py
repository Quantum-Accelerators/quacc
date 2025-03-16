from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from typing import Any


_GEOM_FILE = "geo_end.gen"


def create_dftb_defaults(
    method: Literal["GFN1-xTB", "GFN2-xTB", "DFTB"],
    kpts: tuple | list[tuple] | dict | None = None,
    is_periodic: bool = True,
) -> dict[str, Any]:
    """Create the default calculator kwargs for DFTB+.

    Parameters
    ----------
    method
        Method to use
    kpts
        k-point grid to use
    is_periodic
        Whether the system is periodic

    Returns
    -------
    dict[str, Any]
        Default calculator kwargs
    """
    calc_defaults = {
        "Hamiltonian_": "xTB" if "xtb" in method.lower() else "DFTB",
        "Hamiltonian_MaxSccIterations": 200,
        "kpts": kpts or ((1, 1, 1) if is_periodic else None),
    }
    if "xtb" in method.lower():
        calc_defaults["Hamiltonian_Method"] = method
    return calc_defaults
