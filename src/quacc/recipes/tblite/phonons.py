"""Phonon recipes for TBLite"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import flow
from quacc.recipes.phonons import run_phonons
from quacc.runners.ase import run_calc
from quacc.schemas.ase import summarize_run
from quacc.utils.dicts import merge_dicts

try:
    from tblite.ase import TBLite
except ImportError:
    TBLite = None

if TYPE_CHECKING:
    from typing import Any, Literal

    from ase import Atoms

    from quacc.recipes.phonons import PhononSchema


@flow
def phonon_flow(
    atoms: Atoms,
    method: Literal["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"] = "GFN2-xTB",
    calc_swaps: dict[str, Any] | None = None,
) -> PhononSchema:
    """
    Carry out a single-point calculation.

    Parameters
    ----------
    atoms
        Atoms object
    method
        GFN1-xTB, GFN2-xTB, and IPEA1-xTB.
    calc_swaps
        Dictionary of custom kwargs for the EMT calculator. Set a value to
        `None` to remove a pre-existing key entirely. For a list of available
        keys, refer to the `tblite.ase.TBLite` calculator.

        !!! Info "Calculator defaults"

            ```python
            {"method": method}
            ```
    copy_files
        Files to copy to the runtime directory.

    Returns
    -------
    RunSchema
        Dictionary of results from [quacc.schemas.ase.summarize_run][]
    """

    defaults = {"method": method}
    flags = merge_dicts(defaults, calc_swaps)
    atoms.calc = TBLite(**flags)

    final_atoms = run_phonons(atoms)
    return summarize_run(
        final_atoms,
        input_atoms=atoms,
        additional_fields={"name": "TBLite Static"},
    )
