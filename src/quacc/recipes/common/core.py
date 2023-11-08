"""Core recipes"""
from __future__ import annotations

from typing import TYPE_CHECKING

from quacc import job
from quacc.runners.ase import run_calc

if TYPE_CHECKING:
    from ase.atoms import Atoms
    from ase.calculators.calculator import Calculator
    from numpy.typing import NDArray


@job
def _force_job(atoms: Atoms, calculator: Calculator) -> NDArray:
    """
    Calculate forces

    Parameters
    ----------
    atoms
        Atoms object
    calculator
        Calculator to use

    Returns
    -------
    Forces
        Forces from the calculation.
    """
    atoms.calc = calculator
    return run_calc(atoms).get_forces()
