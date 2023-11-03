"""Core common workflows"""
from typing import TYPE_CHECKING

from quacc import job

if TYPE_CHECKING:
    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from numpy.typing import NDArray

    from quacc.runners.ase import run_calc


@job
def force_job(atoms: Atoms, calculator: Calculator) -> NDArray:
    atoms.calc = calculator
    return run_calc(atoms).get_forces()
